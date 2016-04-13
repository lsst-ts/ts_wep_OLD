from datetime import datetime

import luigi
import numpy as np

from CWFSTask import CWFSTask

class CWFSProcessingTask(luigi.Task):

    date = luigi.DateMinuteParameter(default=datetime.now())

    def output(self):
        return luigi.LocalTarget(
            self.date.strftime("data/cwfsprocessing_%Y%m%d%H%M"))

    def requires(self):
        return CWFSTask()

    def run(self):
        with self.output().open('w') as f:
            f.write(str(1))
        print "Running Processing Task"
