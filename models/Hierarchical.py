import numpy as np
import pandas as pd
import scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import sys

class Hierarchical:

    def __init__(self, data_dir, data_file):
        self.data_dir = data_dir
        self.data_file = data_file
        self.input_path = self.data_dir + '/' + self.data_file
        self.data_matrix = self.create_data_matrix()

    def create_data_matrix(self):
        return pd.read_csv(self.input_path)


def main():
    H = Hierarchical()


if __name__ == "__main__":
    data_dir = sys.argv[1]
    data_file = sys.argv[2]
    main()