import pandas as pd

# edit for snp files with more than one column
def get_snps(filename):
    snp_df = pd.read_csv(filename, squeeze=True)  # squeeze turns it into a pd Series
    return snp_df