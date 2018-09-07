#!/usr/bin/env python

import subprocess
import os
import glob
import sys
import argparse


def main(arguments):
  parser = argparse.ArgumentParser(description=
            '''
            Use system cdo operators to zonally average all CORE forcing with mask to create atmospheric
            forcings for aquaplanet 
            ''',
            epilog = "Written by A.E. Shao 2018")
  parser.add_argument('inpath'  , help='Path to the directory where the CORE fields are stored', type = str)
  parser.add_argument('maskfile', help='Path to file containing the landmask', type = str)
  parser.add_argument('outpath',  help='Path where output will be stored', type = str)
  a = parser.parse_args()

  files = glob.glob(a.inpath + '/*.nc')
  for f in files:
    print("Processing %s" % f)
    subprocess.Popen(['cdo', 'enlarge,'+f, '-zonmean', '-div', f, a.maskfile, a.outpath + '/' + 'aqua_' + os.path.basename(f)]).wait()
if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
