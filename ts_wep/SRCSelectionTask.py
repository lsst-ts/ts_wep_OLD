from datetime import datetime

import luigi

from INSTSignatureRemovalTask import INSTSignatureRemovalTask

class SRCSelectionTask(luigi.Task):

    date = luigi.DateMinuteParameter(default=datetime.now())

    def output(self):
        return luigi.LocalTarget(
            self.date.strftime("data/srcselection_%Y%m%d%H%M"))

    def requires(self):
        return INSTSignatureRemovalTask()

    def run(self):
        with self.output().open('w') as f:
            f.write(str(1))
        print "Running Processing Task"
