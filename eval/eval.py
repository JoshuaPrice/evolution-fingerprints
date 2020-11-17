'''
Usage: python3 eval.py 1.csv 2.csv etc.
The program will compute the evaluation metric for each clustering csv file.
Adjust the following filepath for ground truth labels and raw data if necessary.
'''
LABEL_PATH_INIT = '~/evfi/igsr_samples.tsv'
DATA_PATH_INIT = '/mnt/datadrive/data/full-data-matrix.csv'

import sys

import csv
import numpy as np
import pandas as pd

import pdb

class Eval:

    def __init__(self, label_path, data_path):
        self.label_path = label_path
        self.labels = {}

        self.data_path = data_path
        self.data = None
        self.data_center = None
        self.num_samples = 0
        self.num_features = 0
        	
        self.clusters = []
        self.num_clusters = 0

    def init_label(self):
        '''
        Reads in the ground truth into dict self.labels
        (key) Sample ID > (value) [Population, Superpopulation]
        Population is a subclass of Superpopulation (e.g. FINnish in EURopean)
        '''
        dt_label = pd.read_csv(self.label_path, sep='\t', index_col='Sample name')
        dt_label = dt_label.T.loc[['Population code', 'Superpopulation code']]
        self.labels = dt_label.to_dict('list')

    def init_data(self):
        '''
        Reads in the raw data as a pandas dataframe
        Row index = Sample ID, Column index = SNP ID (rs1234)
        '''
        self.data = pd.read_csv(self.data_path, index_col=0)
        self.num_samples = len(self.data)
        self.num_features = len(self.data.columns)

        dt = self.data.to_numpy()
        self.data_center = np.average(dt, axis=0)

    def read_clusters(self, cluster_path):
        '''
        Reads in the clusters as list of lists
        [[samples in cluster 0], [samples in cluster 1], ..., [samples in cluster k-1]]
        '''
        with open(cluster_path) as cluster_file:
            reader = csv.reader(cluster_file)
            self.clusters = list(reader)
        self.num_clusters = len(self.clusters)

    def calinski_harabasz(self):
        '''
        Calinski-Harabasz score calculates the ratio of:
        (variance between clusters) over (variance within clusters)
        >> [tr(between-matrix) / tr(within-matrix)] * [(num_samples - k) / (k - 1)]
        Bigger C-H score indicates more separation between clusters
        '''
        within_clusters = np.zeros(shape=(self.num_features, self.num_features))
        between_clusters = np.zeros(shape=(self.num_features, self.num_features))

        for cluster_idx in range(self.num_clusters):
            cluster_elems = self.clusters[cluster_idx]
            cluster_data = self.data.loc[cluster_elems].to_numpy()
            cluster_center = np.average(cluster_data, axis=0)
            for elem_idx in range(len(cluster_elems)):
                diff = cluster_data[elem_idx, :] - cluster_center
                within_clusters += np.outer(diff, diff)

            cluster_diff = cluster_center - self.data_center
            between_clusters += len(cluster_elems) * np.outer(cluster_diff, cluster_diff)

            score = (np.trace(between_clusters) / np.trace(within_clusters)) * \
                    ((self.num_samples - self.num_clusters) / (self.num_clusters - 1))

        return np.trace(within_clusters), score

    def get_pop_dist(self, use_superpop=False):
        '''
        Returns a numpy array (num_clusters x num_populations)
        of percentages of each population in each cluster
        '''
        pop_idx_dict = {}
        for pop, superpop in self.labels.values():
            label = superpop if use_superpop else pop
            if label not in pop_idx_dict.keys():
                pop_idx_dict[label] = len(pop_idx_dict.keys())

        pop_dist = np.zeros(shape=(self.num_clusters, len(pop_idx_dict.keys())))
        for cluster_idx in range(self.num_clusters):
            cluster_elems = self.clusters[cluster_idx]
            for elem in cluster_elems:
                pop = self.labels[elem]
                label = pop[1] if use_superpop else pop[0]
                pop_dist[cluster_idx, pop_idx_dict[label]] += 1

        pop_dist = (pop_dist.T / np.sum(pop_dist, axis=1)).T

        return pop_dist


def main(label_path, data_path, cluster_paths):
    ev = Eval(label_path, data_path)
    ev.init_label()
    ev.init_data()
    for cluster_path in cluster_paths:
        ev.read_clusters(cluster_path)
        within_var, ch_score = ev.calinski_harabasz()
        # skip if score is nan (i.e. there is an empty cluster)
        if np.isnan(ch_score):
            continue
        print("Path :", cluster_path, "\nCalinski-Harabasz score", ch_score, "\nwithin_variance", within_var)

        pop_dist = ev.get_pop_dist(use_superpop=True)
        print(pop_dist)

if __name__ == "__main__":
    label_path = LABEL_PATH_INIT
    data_path = DATA_PATH_INIT
    cluster_paths = sys.argv[1:]
    main(label_path, data_path, cluster_paths)

