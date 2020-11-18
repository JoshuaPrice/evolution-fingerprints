import csv
import pandas as pd
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import sys
import utils

class Hierarchical:

    def __init__(self, data_dir, mat_file, set_id, num_clusters):
        self.data_dir = data_dir
        self.data_file = mat_file
        self.set_id = set_id
        self.set_file = set_id + '.csv'
        self.k = num_clusters
        self.input_path = self.data_dir + '/' + self.data_file
        self.snp_set = utils.get_snps(self.data_dir + '/' + self.set_file)
        # self.func_set = utils.get_snps(self.data_dir + '/hap_set.csv')
        self.data_matrix = self.create_data_matrix(self.snp_set)
        self.patients = self.data_matrix.index
        # self.func_data_matrix = self.create_data_matrix(self.func_set)

    def create_data_matrix(self, snp_set):
        full_matrix = pd.read_csv(self.input_path, index_col=0)
        full_mat_cols = full_matrix.columns
        filtered_cols = full_mat_cols.intersection(snp_set)
        # a = full_matrix.loc[:, filtered_cols]
        return full_matrix.loc[:, filtered_cols]

    def cluster_samples(self):
        mat = self.data_matrix
        cluster = AgglomerativeClustering(n_clusters=self.k, affinity='euclidean')
        return cluster.fit_predict(mat)

    def plot_dendogram(self):
        plt.figure(figsize=(10, 7))
        plt.title(self.set_id + " tree")
        dend = hierarchy.dendrogram(hierarchy.linkage(self.data_matrix, method='ward'))
        plt.savefig('../figures/tree.png')

    def save_clusters_to_csv(self):
        cluster_ids = self.cluster_samples()

        # generate empty list of lists
        sample_lists = [[]] * self.k
        for i in range(len(cluster_ids)):
            cluster_id = cluster_ids[i]
            patient_id = self.patients[i]
            sample_lists[cluster_id] = sample_lists[cluster_id] + [patient_id]

        save_path = data_dir + '/clusters_' + self.set_id + '_' + str(self.k) + '.csv'

        with open(save_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(sample_lists)


def main(data_file, mat_file, set_id, num_clusters):
    H = Hierarchical(data_file, mat_file, set_id, num_clusters)
    # clusters = H.cluster_samples()
    # H.plot_dendogram()
    H.save_clusters_to_csv()


if __name__ == "__main__":
    data_dir = sys.argv[1]
    mat_file = sys.argv[2]
    chr_file = sys.argv[3]
    main(data_dir, mat_file, chr_file, 5)