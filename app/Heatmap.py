
from app.settings import *

class Heatmap(object):
    def __init__(self, input, classification, frequency, output, cluster, debug):
        self.input = input
        self.classification = classification
        self.frequency = frequency
        self.output = output
        self.cluster = cluster
        self.debug = debug

        if self.debug:
            logger.setLevel(10)

    def __repr__(self):
        """Returns Heatmap class full object."""
        return "Heatmap({}".format(self.__dict__)

    def run(self):
        # print args
        logger.info(json.dumps(self.__dict__, indent=2))
