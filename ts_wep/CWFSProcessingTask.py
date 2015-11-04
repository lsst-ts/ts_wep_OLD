import luigi
import numpy as np
from CWFSTask import CWFSTask

class CWFSProcessingTask(luigi.Task):

    def output(self):
        return

    def requires(self):
        return CWFSTask()

    def run(self):
        print "Running Processing Task"
