import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

class PlotDist:

    def __init__(self, chr_file):
        self.chr_file = chr_file
        self.chr_data = pd.read_csv(self.chr_file, names=['id','count'])
        self.n = 1233
        self.chr_data['proportion'] = self.chr_data['count'] / self.n

    def plot_and_save(self):
        sns.histplot(data=self.chr_data, x='proportion', bins=20)
        plt.xlabel('Proportion of individuals (n=1233) with a given SNP')
        plt.ylabel('Count')
        plt.title('Histogram of Y-chromosome SNP prevalence in sample')
        plt.savefig('../figures/distplot.png')


if __name__ == "__main__":
    chr_file = sys.argv[1]
    pd = PlotDist(chr_file)
    pd.plot_and_save()