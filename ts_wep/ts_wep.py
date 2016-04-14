# -*- coding: utf-8 -*-

import luigi
from CWFSProcessingTask import CWFSProcessingTask

def main():
    luigi.run(['CWFSProcessingTask', '--workers', '1', '--local-scheduler'])

if __name__ == '__main__':
    main()
