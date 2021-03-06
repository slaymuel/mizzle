""" Print notifications in color

"""

class puts:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    @classmethod
    def warning(cls, string):
        print(cls.FAIL+"Warning"+cls.ENDC+": "+string)
    
    @classmethod
    def success(cls, string):
        print(cls.OKGREEN+"Success"+cls.ENDC+": "+string)