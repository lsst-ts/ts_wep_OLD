import luigi

class INSTSignatureRemovalTask(luigi.Task):

    def output(self):
        return 1

    def requires(self):
        return

    def run(self):
        print "Running Processing Task"
