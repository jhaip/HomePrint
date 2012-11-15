class EndFileException(Exception):
    def __init__(self, args=None):
        self.args = args

class FormatError(Exception):
    def __init__(self, value=None):
        self.value = value
    
    def __str__(self):
        return 'FormatError:' + self.value