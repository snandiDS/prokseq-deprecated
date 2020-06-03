#!/usr/bin/env python
'''----------------------------------------------------------------------------------------
This program is used to convert alignments in BAM file into fastq format
----------------------------------------------------------------------------------------'''

#import built-in modules
import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " RSeQC needs python2.7!\n"
	sys.exit()

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
import sets
from time import strftime

#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyright 2012. All rights reserved."
__credits__ = []
__license__ = "GPL"
__version__="2.6.2"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM format.")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=%default")
	(options,args)=parser.parse_args()

	if not (options.input_file):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file):
		print >>sys.stderr, '\n\n' + input_file + " does NOT exists" + '\n'
		sys.exit(0)

	obj = SAM.ParseBAM(options.input_file)
	obj.stat(q_cut = options.map_qual)


if __name__ == '__main__':
	main()
