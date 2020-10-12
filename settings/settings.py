class LocalSettings:
    """
    Loads settings for local environment, e.g. data_dir (directory containing data).
    """

    def __init__(self):
        self.settings_dict = self.load_local_settings()

    def load_local_settings(self):
        """
        Loads local settings from file and saves them in a dictionary

        :return: dictionary of settings where key = name of setting,
                    value = setting value, e.g {'Name': 'Josh'}
        """
        settings = {}

        settings_file = open("local_settings.txt", "r")
        for line in settings_file:
            key, val = line.split("=")
            settings[key] = val.replace('\n','')
        settings_file.close()

        return settings

    def get_setting(self, key):
        """
        @param key: string representing key of desired setting to retrieve

        :return: value of setting requested. Type can vary.
        """
        return self.settings_dict[key]

