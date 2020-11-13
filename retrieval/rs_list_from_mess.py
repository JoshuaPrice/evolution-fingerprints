import sys
import csv


class FilterRsList:
    '''
    Object to keep track of extracting an rs snp list from a messy file containing substrings separated by
    whitespace. Made this because I copied a table into a messy .txt file.

    Command line usage example:
    python3 rs_list_from_mess.py haps_info_raw.txt hap_set.csv
    '''

    def __init__(self, infoFile):
        self.filename = infoFile
        self.rs_list = []

    def extractString(self):
        '''
        Extracts a single string with all contents of any file.
        :return: Single string with all contents of the file.
        '''
        contents_str = ''
        with open(self.filename, "r") as f:
            for line in f:
                contents_str += line

        return contents_str

    def parseRs(self):
        '''
        Parses the already-extracted string into a list consisting only of substrings starting with rs.
        :return: List of rs snp names.
        '''
        contents_str = self.extractString()
        str_tokens = contents_str.split()
        self.rs_list = [token for token in str_tokens if token[:2]=='rs']
        return self.rs_list

    def saveToCsv(self, outpath):
        '''
        Saves already-computed rs snp list to the designated file.
        :param outpath: Filepath to save snp list .csv to.
        :return: None
        '''
        with open(outpath, 'w') as f:
            writer = csv.writer(f)
            for rs in self.rs_list:
                writer.writerow([rs])


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_path = sys.argv[2]

    filter_obj = FilterRsList(input_file)
    rs_list = filter_obj.parseRs()
    print(rs_list)

    filter_obj.saveToCsv(output_path)


