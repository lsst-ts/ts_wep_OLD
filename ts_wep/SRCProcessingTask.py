import luigi
from SRCSelectionTask import SRCSelectionTask

class SRCProcessingTask(luigi.Task):

    def output(self):
        return 1

    def requires(self):
        return SRCSelectionTask()

    def run(self):
        print "Running Processing Task"
