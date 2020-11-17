import sys
import pandas as pd


class MergeMatrices:

    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.chr_list = ['Y','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
        self.full_data_matrix = self.merge_all()

    def merge_all(self):
        first_chr = self.chr_list[0]
        chr_filename = self.data_dir + '/' + 'matrix-chr' + first_chr + '.csv'
        running_mat = pd.read_csv(chr_filename, index_col=0)

        for i in range(1,len(self.chr_list)):
            chr = self.chr_list[i]
            chr_filename = self.data_dir + '/' + 'matrix-chr' + chr + '.csv'
            add_mat = pd.read_csv(chr_filename, index_col=0)
            running_mat = running_mat.merge(add_mat, left_index=True, right_index=True)

        return running_mat

    def save_full_mat(self):
        self.full_data_matrix.to_csv(self.data_dir + "/full-data-matrix.csv")


if __name__ == '__main__':
    data_dir = sys.argv[1]
    mm = MergeMatrices(data_dir)
    mm.save_full_mat()
