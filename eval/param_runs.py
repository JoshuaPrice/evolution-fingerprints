from eval import *
import pandas as pd
import matplotlib.pyplot as plt
import pdb

KS = [2, 3, 4, 5, 6]
PS = [2, 3, 5, 10]
RUNS = [0, 1, 2]

ev = Eval('/mnt/datadrive/data/full-data-matrix.csv')
'''
for k in KS:
    for p in PS:
        if k == 2 and p == 3:
            continue
        print('k =', k, 'p =', p, end='\t')
        score_list = []
        for run in RUNS:
            filename = '../smalldata/kmeans_clusters/clusters_kmeans_hap_set_k' + str(k) + \
                    '_p' + str(p) + '_' + str(run) + '.csv'
            ev.read_clusters(filename)
            silhouette_score = ev.silhouette()
            print('Sil', format(silhouette_score, ".3f"), end=' ')

            within_var, ch_score = ev.calinski_harabasz()
            print('CH', format(ch_score, ".2f"), 'WCSS', format(within_var, ".2f"), end='\t')
            # soscore = ev.sum_of_squares()
            # score_list.append(soscore)
        # print(format(min(score_list), ".2f"), end='\t')
        print()
'''
file_pop = '~/evfi/igsr_samples.tsv'
#filedir = '../smalldata/clusters_'
filedir = '../smalldata/kmeans_clusters/final_clusters/best_clusters_kmeans_'
filenames = ['ADHD_k3_p16_1.csv', 'Alzheimers_k3_p8_2.csv', 'Coronary_Artery_k3_p24_2.csv',\
        'Hapset_k3_p5_1.csv', 'Height_k3_p40_1.csv']
#filenames = ['ADHD_set_3.csv', 'alzheimers_set_3.csv', 'coronary_artery_3.csv',\
#        'hap_set_3.csv', 'height_set_3.csv']
for filename in filenames:
    filepath = filedir + filename
    ev.read_clusters(filepath)
    ev.init_pop(file_pop)
    pop_dist = ev.get_pop_dist()

    plt.figure(figsize=(3,2))
    pop_dist.plot(kind = 'bar')

    # label axes and title
    plt.xlabel('Cluster')
    plt.ylabel('Percentage of population in cluster')
    if filename == 'Hapset_k3_p5_1.csv':
        cat = 'haplotype_set'
    elif filename == 'Coronary_Artery_k3_p24_2.csv':
        cat = 'coronary_artery_disease'
    else:
        cat = filename.split('_')[0]
    plt.title('Population distribution for ' + cat + ' (K-means)') 
    
    plt.tight_layout()
    plt.savefig('../figures/' + cat + 'popdist_kmeans.png')




