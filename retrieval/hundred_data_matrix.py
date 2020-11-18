import pandas as pd
import random
import sys


class SampleDataMat:

    def __init__(self, in_file):
        self.in_file = in_file
        self.small_mat = self.make_small_mat()

    def make_small_mat(self):
        full_mat = pd.read_csv(self.in_file)
        n =full_mat.shape[0]
        random_indices = random.sample(range(n), 20)
        hundred_mat = full_mat.iloc[random_indices,:]
        return hundred_mat.set_index('Unnamed: 0')

    def save_small_to(self, out_file):
        self.small_mat.to_csv(out_file)


if __name__ == '__main__':
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    sdm = SampleDataMat(in_file)
    sdm.save_small_to(out_file)