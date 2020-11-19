import numpy as np
import pandas as pd
import random
import csv
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import utils


random.seed()

# change this to your data file path
all_matrix_filepath = "smalldata/full-data-matrix.csv"

hapset_filepath = "smalldata/hap_set.csv"
ADHD_filepath = "smalldata/ADHD.csv"
alzheimers_filepath = "smalldata/alzheimers_set.csv"
artery_filepath = "smalldata/coronary_artery.csv"
height_filepath = "smalldata/height.csv"

patient_region_filepath = "/Users/juliapark/GitHub/julia-evolution-data/igsr_samples.tsv"

chosen_k = 4
hapset_p = [5, 3, 7] # by scree plot, by .8, by .9
adhd_p = [7, 16, 21]
alzheimers_p = [8, 8, 10]
artery = [3, 24, 30]
height = [5, 40, 51]




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

    def __init__(self, k=5, predictor_count=None, max_iter=10000, tol=1e-5):
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


    def read_data_frequent(self, csv_path):
        """
        OLD MODEL: DO NOT USE
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


        # by frequency
        frequency = np.sum(inputs, axis=0)
        ind = np.argpartition(frequency, -1 * self.predictor_count)[-1 * self.predictor_count:]
        ind = ind[np.argsort(frequency[ind])]

        final_features = inputs[:, ind]
        feature_names = feature_names[ind]

        return patient_ID, final_features, feature_names


    def fit(self, data, patient_ID):
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

def create_data_matrix(snp_set):
    full_matrix = pd.read_csv(all_matrix_filepath, index_col=0)
    full_mat_cols = full_matrix.columns
    filtered_cols = full_mat_cols.intersection(snp_set)
    data_matrix = full_matrix.loc[:, filtered_cols]
    data_matrix = data_matrix.to_numpy()

    patient_ID = np.loadtxt(all_matrix_filepath, delimiter=',', dtype="U10", skiprows=1, usecols=0)
    return data_matrix, patient_ID


def transform_data_PCA(data_matrix, n_components):
    pca = PCA(n_components=n_components)
    final_features = pca.fit_transform(data_matrix)
    return final_features


def run_experiment(filter_filepath, set_type, ks, predictor_counts, rounds=3):
    # snp_set = utils.get_snps(filter_filepath)
    snp_set = np.loadtxt(filter_filepath, dtype="U10", delimiter=',', usecols=1)
    data_matrix, patient_ID = create_data_matrix(snp_set)
    for k in ks:
        for predictor_count in predictor_counts:
            features = transform_data_PCA(data_matrix, predictor_count)
            for i in range(rounds):
                model = Kmeans(k=k)

                model.fit(features, patient_ID)

                output_filepath = "clusters_kmeans_%s_k%s_p%s_%s.csv" % (set_type, k, predictor_count, i)
                w = csv.writer(open(output_filepath, "w"))
                for key, val in model.patient_clusters.items():
                    w.writerow(val)
                print("finished ", output_filepath)


def find_best_predictor_count(inputs, set_type):
    pca = PCA()
    features = pca.fit_transform(inputs)
    i = 0
    explained = 0
    for p in pca.explained_variance_ratio_:
        explained += p
        i += 1
        if explained > .75:
            break
    output = "%s components needed for .8 variance" % i
    print(output)

    PC_values = np.arange(pca.n_components_) + 1
    plt.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=1)
    title = '%s Scree Plot' % (set_type)
    plt.title(title)
    plt.xlabel('component')
    plt.ylabel('explained variance')
    plt.show()



def main():
    all_sets_filepath = [hapset_filepath, ADHD_filepath, alzheimers_filepath, artery_filepath, height_filepath]
    all_sets_types = ["Hapset", "ADHD", "Alzheimers", "Coronary_Artery", "Height"]
    # uncomment this block to get scree plot
    '''
    for i in range(0, 5):
        set_type = all_sets_types[i]
        snp_set = np.loadtxt(all_sets_filepath[i], dtype="U10", delimiter=',', usecols=1)
        data, patient_ID = create_data_matrix(snp_set)
        find_best_predictor_count(data, set_type)
    '''

    # '''
    # uncomment this block to cluster
    ks = [3]
    all_sets_predictor = [[5], [16], [8], [24], [40]]
    for i in range(0, 5):
        run_experiment(all_sets_filepath[i], all_sets_types[i], ks, all_sets_predictor[i])
    # '''

    df = pd.read_csv(patient_region_filepath, sep="\t")
    # df = df[['Sample name', 'Sex', 'Population code', 'Population name', 'Superpopulation code', 'Superpopulation name']]
    # cluster0 = df[df['Sample name'].isin(model.patient_clusters[0])]

    return







if __name__ == '__main__':
    main()