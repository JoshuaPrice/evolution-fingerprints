import numpy as np
import pandas as pd
import random
import csv
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import utils


random.seed()

# change this to your data file path
chrY_matrix_filepath = "/Users/juliapark/GitHub/julia-evolution-data/matrix-chrY.csv"
patient_region_filepath = "/Users/juliapark/GitHub/julia-evolution-data/igsr_samples.tsv"
all_matrix_filepath = "/Users/juliapark/GitHub/julia-evolution-data/full-data-matrix.csv"
hapset_filepath = "/Users/juliapark/GitHub/evolution-fingerprints/smalldata/hap_set.csv"

chosen_k = 4
hapset_p = 5


class Kmeans:
    """Kmeans clustering.

    Properties:
         k: the number of clusters. (int)
            predictor_count: number of predictors to use (int)
            max_iter: maximum number of iterations (int)
            tol: used for testing convergence (float)
            centroids: values of centroids (dictionary with keys 0 to k-1, values arrays of size predictor_count)
            clusters: the data points assigned to each cluster (dictionary with keys 0 to k-1, values arrays of predictor arrays)
            patient_clusters: the patients assigned to each cluster (dictionary with keys 0 to k-1, values arrays of patient IDs)
    """

    def __init__(self, k=5, predictor_count=20, max_iter=10000, tol=1e-5):
        """
        Args:
            k: the number of clusters. (int)
            predictor_count: number of predictors to use (int)
            max_iter: maximum number of iterations (int)
            tol: used for testing convergence (float)
        """
        self.k = k
        self.predictor_count = predictor_count
        self.max_iter = max_iter
        self.tol = tol
        self.centroids = {}
        self.clusters = {}
        self.patient_clusters = {}

    def read_data(self, csv_path):
        """
        Args:
            csv_path: path to a csv file
        """
        # Load features
        headers = np.genfromtxt(csv_path, delimiter=',', dtype="U10", max_rows=1)

        x_cols = np.asarray((np.arange(1, len(headers))))
        patient_cols = [0]

        inputs = np.loadtxt(csv_path, delimiter=',', skiprows=1, usecols=x_cols)
        patient_ID = np.loadtxt(csv_path, delimiter=',', dtype="U10", skiprows=1, usecols=patient_cols)
        feature_names = headers[1:]

        # by PCA
        pca = PCA(n_components=self.predictor_count)
        final_features = pca.fit_transform(inputs)
        feature_names = pca

        # by frequency
        #frequency = np.sum(inputs, axis=0)
        #ind = np.argpartition(frequency, -1 * self.predictor_count)[-1 * self.predictor_count:]
        #ind = ind[np.argsort(frequency[ind])]

        #final_features = inputs[:, ind]
        #feature_names = feature_names[ind]

        return patient_ID, final_features, feature_names


    def fit(self, patient_ID, data):
        """
        Args:
            patient_ID: vector of patient IDs
            data: snp data matrix
        """
        # random initialization
        for i in range(self.k):
            self.centroids[i] = data[random.randrange(start=0, stop=len(patient_ID)-1)]

        # static initialization
        #for i in range(self.k):
        #   self.centroids[i] = data[i]

        for i in range(self.max_iter):
            for j in range(self.k):
                self.clusters[j] = []
                self.patient_clusters[j] = []

            for k in range(data.shape[0]):
                x = data[k]
                distances = [np.linalg.norm(x - self.centroids[centroid]) for centroid in self.centroids]
                closest_class = distances.index(min(distances))
                self.clusters[closest_class].append(x)
                self.patient_clusters[closest_class].append(patient_ID[k])

            prev_centroids = self.centroids

            optimized = True

            for cluster in self.clusters:
                self.centroids[cluster] = np.average(self.clusters[cluster], axis=0)
                if np.sum((self.centroids[cluster] - prev_centroids[cluster]) / self.centroids[cluster] * 100.0) > self.tol:
                    optimized = False

            if optimized:
                break

def read_inputs(data_filepath, snp_set):
    full_matrix = pd.read_csv(data_filepath, index_col=0)

    full_mat_cols = full_matrix.columns
    filtered_cols = full_mat_cols.intersection(snp_set)
    patient_ID = np.loadtxt(data_filepath, delimiter=',', dtype="U10", skiprows=1, usecols=0)
    return full_matrix.loc[:, filtered_cols], patient_ID

def run_experiment(data_filepath, filter_filepath, ks, predictor_counts, rounds=3):
    snp_set = utils.get_snps(filter_filepath)
    features, patient_ID = read_inputs(data_filepath, snp_set)
    features = features.to_numpy()
    for k in ks:
        for predictor_count in predictor_counts:
            for i in range(rounds):
                model = Kmeans(k=k, predictor_count=predictor_count)
                model.fit(patient_ID, features)

                output_filepath = "clusters_kmeans_hap_set_k%s_p%s_%s.csv" % (k, predictor_count, i)
                w = csv.writer(open(output_filepath, "w"))
                for key, val in model.patient_clusters.items():
                    w.writerow(val)
                print("finished ", output_filepath)


def find_best_predictor_count(inputs, set_type):
    pca = PCA(n_components=50)
    features = pca.fit_transform(inputs)
    PC_values = np.arange(pca.n_components_) + 1
    plt.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=1)
    title = '%s Scree Plot' % (set_type)
    plt.title(title)
    plt.xlabel('component')
    plt.ylabel('explained variance')
    plt.show()


def main():
    # data = read_inputs(all_matrix_filepath, hapset_filepath)
    # find_best_predictor_count(data, )
    ks = [2, 3, 4, 5, 6]
    predictor_counts = [2, 3, 5, 10]
    run_experiment(all_matrix_filepath, hapset_filepath, ks, predictor_counts)


    # df = pd.read_csv(patient_region_filepath, sep="\t")
    # df = df[['Sample name', 'Sex', 'Population code', 'Population name', 'Superpopulation code', 'Superpopulation name']]
    # cluster0 = df[df['Sample name'].isin(model.patient_clusters[0])]

    return







if __name__ == '__main__':
    main()