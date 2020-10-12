class LocalSettings:
    """
    Loads settings for local environment, e.g. data_dir (directory containing data).
    Example usage:
        import settings
        s = settings.LocalSettings()
        data_dir = s.get_setting('data_dir')
    """

    def __init__(self, filename='local_settings.txt'):
        self.filename=filename
        self.settings_dict = self.load_local_settings()

    def load_local_settings(self):
        """
        Loads local settings from file and saves them in a dictionary

        :return: dictionary of settings where key = name of setting,
                    value = setting value, e.g {'Name': 'Josh'}
        """
        settings = {}

        settings_file = open(self.filename, "r")
        for line in settings_file:
            if line not in ['\n', '\r\n']:
                k, v = line.split("=")
                settings[k] = v.replace('\n','')
        settings_file.close()

        return settings

    def get_setting(self, key):
        """
        @param key: string representing key of desired setting to retrieve

        :return: value of setting requested. Type can vary.
        """
        return self.settings_dict[key]


if __name__ == '__main__':
    # if called at command line, print all local settings
    s = LocalSettings()
    print("Local settings are:\n")
    for key in s.settings_dict:
        print(key + ": " + s.get_setting(key))
    print("\n")
