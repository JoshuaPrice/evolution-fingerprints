import pandas as pd

# edit for snp files with more than one column
def get_snps(filename):
    snp_df = pd.read_csv(filename, names=['chr','id'])  # squeeze turns it into a pd Series
    return snp_df['id']
