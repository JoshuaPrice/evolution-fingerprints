# from evfi.Settings.settings import LocalSettings

class SimData:

    def __init__(self, directory, num_haplotypes=5, num_markers=3):
        self.directory = directory
        self.num_haplotypes = num_haplotypes
        self.num_markers = num_markers