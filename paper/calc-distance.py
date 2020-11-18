import pandas as pd
import sys
import csv

class CalcDist:

    def __init__(self, cluster1_file, cluster2_file):
        self.cluster1_file = cluster1_file
        self.cluster2_file = cluster2_file
        self.cluster1_df = self.file_to_df(cluster1_file)
        self.cluster2_df = self.file_to_df(cluster2_file)

    def file_to_df(self, cluster_file):
        cluster_df = pd.DataFrame(columns=['id','cluster'])
        with open(cluster_file) as f:
            cluster_rows = csv.reader(f, delimiter=',')
            cluster_num = 1
            for row in cluster_rows:
                for id in row:
                    cluster_df = cluster_df.append({'id': id, 'cluster': cluster_num}, ignore_index=True)
                cluster_num += 1

        return cluster_df

    def get_distance(self):
        running_dists = []
        for idx, row in self.cluster1_df.iterrows():
            cluster1_id = row['cluster']
            cluster1_bool = self.cluster1_df['cluster'] == cluster1_id

            cluster2_id = row['cluster']
            cluster2_bool = self.cluster2_df['cluster'] == cluster2_id

            matching_bool = cluster1_bool == cluster2_bool
            matching_ratio = sum(matching_bool) / len(matching_bool)
            running_dists.append(matching_ratio)

        overall_dist = 1 - sum(running_dists) / len(running_dists)

        return overall_dist


if __name__ == '__main__':
    # example usage: python3 calc-distance.py clusters_hap_set_3.csv clusters_ADHD_set_3.csv
    cluster1_file = sys.argv[1]
    cluster2_file = sys.argv[2]
    cd = CalcDist(cluster1_file, cluster2_file)
    dist = cd.get_distance()
    print("overall distance is " + str(dist))
