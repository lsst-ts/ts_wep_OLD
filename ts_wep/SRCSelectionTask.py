import luigi
from INSTSignatureRemovalTask import INSTSignatureRemovalTask

class SRCSelectionTask(luigi.Task):

    def output(self):
        return 1

    def requires(self):
        return INSTSignatureRemovalTask()

    def run(self):
        print "Running Processing Task"
