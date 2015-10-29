# -*- coding: utf-8 -*-

import luigi
from CWFSTask import CWFSTask

def main():
    luigi.run(['CWFSTask', '--workers', '1', '--local-scheduler'])

if __name__ == '__main__':
    main()
