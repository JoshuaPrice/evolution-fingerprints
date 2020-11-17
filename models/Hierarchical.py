import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import sys
import utils

class Hierarchical:

    def __init__(self, data_dir, data_file):
        self.data_dir = data_dir
        self.data_file = data_file
        self.input_path = self.data_dir + '/' + self.data_file
        self.hap_set = utils.get_snps(self.data_dir + '/hap_set.csv')
        self.func_set = utils.get_snps(self.data_dir + '/hap_set.csv')
        self.evol_data_matrix = self.create_data_matrix(self.hap_set)
        # self.func_data_matrix = self.create_data_matrix(self.func_set)

    def create_data_matrix(self, ):
        full_matrix = pd.read_csv(self.input_path, index_col=0)
        full_mat_cols = full_matrix.columns
        filtered_cols = full_mat_cols.intersection(self.hap_set)
        a = full_matrix.loc[:, filtered_cols]
        return full_matrix.loc[:, filtered_cols]

    def cluster_samples(self, mat):
        cluster = AgglomerativeClustering(n_clusters=3, affinity='euclidean')
        return cluster.fit_predict(mat)

    def plot_evol_dendogram(self):
        plt.figure(figsize=(10, 7))
        plt.title("Evoluationary (null hypothesis) tree")
        dend = hierarchy.dendrogram(hierarchy.linkage(self.data_matrix, method='ward'))
        plt.savefig('figures/evolutionary_tree.png')

    def plot_func_dendogram(self):
        pass

def main(data_file, mat_file):
    H = Hierarchical(data_file, mat_file)
    H.plot_evol_dendogram()


if __name__ == "__main__":
    data_dir = sys.argv[1]
    mat_file = sys.argv[2]
    main(data_dir, mat_file)