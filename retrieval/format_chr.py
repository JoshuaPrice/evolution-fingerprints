import pandas as pd
import sys
import os


class FormatChr:
    '''
    Object to read in data directly from a 1000 Genomes file and save a formatted, filtered, and transposed data matrix.
    example command line usage: python3 format_chr.py "/shared/data" "chrY.csv"
    In other programs (in this same directory), can directly set up pandas dataframes by calling:
        import format_chr
        format_obj = format_chr.FormatChr(data_dir, data_file)
        format_obj.write_readable_file()
        data_matrix = format_obj.read_and_format
    '''

    def __init__(self, data_dir, chr_id):
        self.input_file = 'chr' + str(chr_id) + '.csv'
        self.data_dir = data_dir
        self.temp_readable_file = data_dir + '/temp_' + self.input_file
        self.interesting_snps = self.get_interesting_snps(chr_id)

    def get_interesting_snps(self, chr_id):
        chr_file = self.data_dir + '/all_sets.csv'
        all_interesting_snps_df = pd.read_csv(chr_file, names=['chr', 'id'], index_col=None)
        chr_intereting_snps = all_interesting_snps_df[all_interesting_snps_df['chr']==chr_id]
        return chr_intereting_snps['id'].to_list()

    def write_readable_file(self):
        '''
        Pandas can only read in consistently-dimensioned files. Needs the same number of tabs in each row.
        This method removes all lines without the correct number of tabs.
        :return: None
        '''
        raw_data_file = open(self.data_dir + '/' + self.input_file, 'r')
        temp_readable_file = open(self.temp_readable_file, 'w')
        read_in_line = False
        tab_count_constraint = 0
        for line in raw_data_file:
            if line[:6] == "#CHROM":
                read_in_line = True
                tab_count_constraint = line.count('\t')
            line_filtered = line.replace(".", "0.0").replace("0|0", "0").replace("0|1", "0.5").replace("1|0", "0.5").replace("1|1", "1.0")

            consider_line_if_includes = self.interesting_snps + ["ID"] # need to include header line
            if read_in_line and line.count('\t') == tab_count_constraint and any(interest in line for interest in consider_line_if_includes):
                temp_readable_file.write(line_filtered)

        raw_data_file.close()
        temp_readable_file.close()

    def read_and_format(self):
        '''
        Reads the data file into pandas, removes the columns we don't want, replaces all strange values with 0,
        filters to only include named SNPs (rs...), transposes matrix so subjects are rows and SNPs are columns,
        saves filtered & transposed matrix to

        :return: pandas dataframe of 0s/1s with rows being samples (people), cols being SNP IDs
        '''
        # read and delete the temp file we made
        mat_original = pd.read_csv(self.temp_readable_file, sep='\t')
        os.remove(self.temp_readable_file)

        # filter to only include SNPs with rs IDs and remove non-binary columns
        labeled_snps = mat_original[mat_original["ID"] != "0"]
        labeled_binary = labeled_snps.drop(columns=["#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])

        # set any weird values in the matrix to 0, make all values ints (0 or 1) instead of strings
        for col in labeled_binary.columns[1:]:
            labeled_binary.loc[~labeled_binary[col].isin([0, 0.0, 0.5, 1.0, 1]), col] = 0.0
        labeled_binary = labeled_binary.set_index("ID").astype(float)

        # transpose to make Julia happy
        labeled_flipped = labeled_binary.transpose()
        return labeled_flipped

    def format_and_save(self):
        '''
        Saves file containing matrix data as well as csv containing SNPs ordered by frequency.
        :return: None
        '''
        labeled_flipped = self.read_and_format()
        print("got past read_and_format")
        labeled_flipped.to_csv(self.data_dir + "/matrix-" + self.input_file)

        labeled_binary = labeled_flipped.transpose()
        labeled_snp_counts = labeled_binary.sum(axis=1)
        sorted_idx = labeled_snp_counts.sort_values(ascending=False)
        sorted_idx.to_csv(self.data_dir + "/sorted-snps-" + self.input_file)


def main(data_dir, data_file):
    format_obj = FormatChr(data_dir, data_file)
    print("Created FormatChr object. Now writing temp file")
    format_obj.write_readable_file()
    print("Finished writign temp file. Now formatting and saving.")
    format_obj.format_and_save()
    print("Exiting program.")


if __name__ == "__main__":
    # example usage: python3 format_chr.py "/shared/data" "Y"
    data_dir = sys.argv[1]
    chr_id = sys.argv[2]
    main(data_dir, chr_id)
