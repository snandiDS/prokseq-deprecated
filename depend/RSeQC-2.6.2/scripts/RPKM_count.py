#!/usr/bin/env python

'''--------------------------------------------------------------------------------
calculate raw read, RPKM for each exon, intron, mRNA regions in provided gene list
--------------------------------------------------------------------------------'''

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
from qcmodule import bam_cigar
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


def printlog (mesg):
	'''print progress into stderr and log file'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	LOG=open('class.log','a')
	print >>sys.stderr,mesg
	print >>LOG,mesg

def load_chromsize(file):
        '''read chrom.size file'''
        chromSize={}
        for line in open(file,'r'):
                if line.startswith('#'):continue
                if not line.strip():continue
                fields = line.strip().split()
                chromSize[fields[0]] = int(fields[1])
        return chromSize

def build_range(refgene):
	'''build ranges for exonic region'''
	ranges={}
	for line in open(refgene,'r'):
		try:
			if line.startswith(('#','track','browser')):continue

			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0].upper()
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5].replace(" ","_")
				
			exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
		except:
			print >>sys.stderr,"[NOTE:input bed must be 12-column] skipped this line: " + line,
			continue	

		for st,end in zip(exon_starts,exon_ends):	
			if chrom not in ranges:
				ranges[chrom] = Intersecter()
			ranges[chrom].add_interval( Interval( st, end ) )
	return ranges
	

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM format (SAM is not supported). [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed fomat. [required]")
	parser.add_option("-d","--strand",action="store",type="string",dest="strand_rule",default=None,help="How read(s) were stranded during sequencing. For example: --strand='1++,1--,2+-,2-+' means that this is a pair-end, strand-specific RNA-seq, and the strand rule is: read1 mapped to '+' => parental gene on '+'; read1 mapped to '-' => parental gene on '-'; read2 mapped to '+' => parental gene on '-'; read2 mapped to '-' => parental gene on '+'.  If you are not sure about the strand rule, run \'infer_experiment.py' default=%default (Not a strand specific RNA-seq data)")
	parser.add_option("-u","--skip-multi-hits",action="store_true",dest="skip_multi",help="How to deal with multiple hit reads. Presence this option renders program to skip multiple hits reads.")
	parser.add_option("-e","--only-exonic",action="store_true",dest="only_exon",help="How to count total reads. Presence of this option renders program only used exonic (UTR exons and CDS exons) reads, otherwise use all reads.")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be called \"uniquely mapped\". default=%default")
	
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file and options.refgene_bed):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file + '.bai'):
		print >>sys.stderr, "cannot find index file of input BAM file"
		print >>sys.stderr, options.input_file + '.bai' + " does not exists"
		sys.exit(0)
	for file in (options.input_file, options.refgene_bed):
		if not os.path.exists(file):
			print >>sys.stderr, file + " does NOT exists" + '\n'
			sys.exit(0)

	obj = SAM.ParseBAM(options.input_file)
	OUT = open(options.output_prefix + '_read_count.xls','w')

	#++++++++++++++++++++++++++++++++++++determine strand rule
	strandRule={}
	if options.strand_rule is None:													# Not strand-specific
		pass																
	elif len(options.strand_rule.split(',')) ==4:									#PairEnd, strand-specific
		for i in options.strand_rule.split(','):strandRule[i[0]+i[1]]=i[2]
	elif len(options.strand_rule.split(',')) ==2:									#singeEnd, strand-specific
		for i in options.strand_rule.split(','):strandRule[i[0]]=i[1]
	else:
		print >>sys.stderr, "Unknown value of option :'strand_rule' " +  options.strand_rule
		sys.exit(1)	

	#++++++++++++++++++++++++++++++++++++counting reads
	print >>sys.stderr, "Retrieve exon regions from  "+ options.refgene_bed + '...'
	gene_ranges = build_range( options.refgene_bed)
	#print gene_ranges['ERCC-00002'].find(0,100)
	print >>sys.stderr, "Counting total reads ... ",
	
	total_reads =0
	total_tags =0
	total_exonic_tags =0
	
	try:
		while(1):
			aligned_read = obj.samfile.next()
			if aligned_read.is_qcfail:continue			#skip low quanlity					
			if aligned_read.is_duplicate:continue		#skip duplicate read
			if aligned_read.is_secondary:continue		#skip non primary hit
			if aligned_read.is_unmapped:continue		#skip unmap read
			if options.skip_multi:
				if aligned_read.mapq < options.map_qual:
					continue
			total_reads +=1
			chrom = obj.samfile.getrname(aligned_read.tid).upper()
			hit_st = aligned_read.pos
			exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)	
			total_tags += len(exon_blocks)

			for exn in exon_blocks:	
				mid = exn[1] + int((exn[2]-exn[1])/2)
				#print chrom,mid,mid+1
				#print gene_ranges[chrom].find(mid,mid+1)
				if (chrom in gene_ranges) and len(gene_ranges[chrom].find(mid,mid+1)) >0:
					total_exonic_tags += 1
					
	except StopIteration:
		print >>sys.stderr, "Done"
	print >>sys.stderr, "Total Reads = %-20s" % (str(total_reads))
	print >>sys.stderr, "Total Tags = %-20s" % (str(total_tags))
	print >>sys.stderr, "Total Exon Tags = %-20s" % (str(total_exonic_tags))

	if total_tags >0 and total_exonic_tags>0:
		if options.only_exon:
			denominator = total_exonic_tags
		else:
			denominator = total_tags
	else:
		print >>sys.stderr, "Total tags cannot be 0 or negative number"
		sys.exit(1)
	
	#++++++++++++++++++++++++++++++++++++++++++++++++
	obj = SAM.ParseBAM(options.input_file)
	if options.strand_rule is None:
		OUT.write('#chrom' + '\t' + 'st' + '\t' + 'end' + '\t' + 'accession' + '\t' + 'score' + '\t' + 'gene_strand' + '\t' + 'tag_count' + '\t' + 'RPKM' + '\n')
	else:
		OUT.write('#chrom' + '\t' + 'st' + '\t' + 'end' + '\t' + 'accession' + '\t' + 'score' + '\t' + 'gene_strand' + '\t' + 'tag_count_Forward' + '\t' + 'tag_count_Reverse' +'\t' + 'RPKM_Forward' + '\t' + 'RPKM_Reverse' + '\n')
	genome_total_read=0
	genome_unique_read=0
	gene_finished=0
	#calculate raw count, RPKM for each gene
	for line in open(options.refgene_bed,'r'):
		exon_range=Intersecter()
		intron_range=Intersecter()		
		if line.startswith(('#','track','browser')):continue   
		fields = line.split()
		chrom     = fields[0]
		tx_start  = int( fields[1] )
		tx_end    = int( fields[2] )
		geneName      = fields[3]
		gstrand    = fields[5].replace(" ","_")
		cds_start = int( fields[6] )
		cds_end   = int( fields[7] )
	    	
		exon_starts = map( int, fields[11].rstrip( ',\n' ).split( ',' ) )
		exon_starts = map((lambda x: x + tx_start ), exon_starts)
		exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
		exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
		intron_starts = exon_ends[:-1]
		intron_ends = exon_starts[1:]
		
		plus_ranges=Intersecter()
		minus_ranges=Intersecter()
		unstrand_ranges= Intersecter()
		
		try:
			alignedReads = obj.samfile.fetch(chrom,tx_start,tx_end)
		except:
			print >>sys.stderr, "No alignments for " + geneName + ". Skip"
			continue
		for aligned_read in alignedReads:
			flag=0
			if aligned_read.is_qcfail:continue			#skip low quanlity					
			if aligned_read.is_duplicate:continue		#skip duplicate read
			if aligned_read.is_secondary:continue		#skip non primary hit
			if aligned_read.is_unmapped:continue		#skip unmap read
			
			if options.skip_multi:
				if aligned_read.mapq < options.map_qual:
					continue

			if aligned_read.is_paired:						#pair end
				if aligned_read.is_read1:read_id = '1'
				if aligned_read.is_read2:read_id = '2'
			else:read_id = ''								#single end
			
			if aligned_read.is_reverse:map_strand = '-'
			else:map_strand = '+'				
			strand_key = read_id + map_strand				#used to determine if a read should assign to gene(+) or gene(-)

			hit_st = aligned_read.pos
			exon_blocks = bam_cigar.fetch_exon(chrom, hit_st, aligned_read.cigar)
						
			#construct bitset
			if options.strand_rule is not None:	
				if strandRule[strand_key] == '+':
					for block in exon_blocks:
						mid = block[1] + int((block[2] - block[1])/2)
						plus_ranges.add_interval( Interval( mid,mid+1 ) )
				elif strandRule[strand_key] == '-':
					for block in exon_blocks:
						mid = block[1] + int((block[2] - block[1])/2)
						minus_ranges.add_interval( Interval( mid,mid+1 ) )
			elif options.strand_rule is None:	
				for block in exon_blocks:
					mid = block[1] + int((block[2] - block[1])/2)
					unstrand_ranges.add_interval( Interval( mid,mid+1 ) )
		mRNA_plus_hits =0
		mRNA_plus_rpkm =0.0
		
		mRNA_minus_hits =0
		mRNA_minus_rpkm =0.0
		
		mRNA_hits =0
		mRNA_rpkm =0.0
		
		mRNA_length=0
		
		#assign reads to region:exon,intron,mRNA
		if (options.strand_rule is not None):	#this is strand specific
			if gstrand == '-':
				intronNum=len(intron_starts)	
				exonNum=len(exon_starts)
			elif gstrand == '+':
				intronNum=1
				exonNum=1			
			#assign reads to intron regions
			for st,end in zip(intron_starts,intron_ends):
				if end >st:
					size = end - st
				elif end == st:
					size = 1
				hits_plus = len(plus_ranges.find(st,end))
				hits_minus = len(minus_ranges.find(st,end))
				hits_plus_rpkm = hits_plus*1000000000.0/(size*denominator)
				hits_minus_rpkm = hits_minus*1000000000.0/(size*denominator)
				print >>OUT, '\t'.join(['%s','%d','%d','%s','%d','%s','%d','%d','%.3f','%.3f']) % (chrom,st,end,geneName + "_intron_" + str(intronNum),0,gstrand,hits_plus,hits_minus,hits_plus_rpkm,hits_minus_rpkm)
				if gstrand == '-':intronNum -= 1
				elif gstrand == '+':intronNum +=1
			#assign reads to exon regions
			for st,end in zip(exon_starts,exon_ends):
				if end >st:
					size = end - st
				elif end == st:
					size = 1			
				hits_plus = len(plus_ranges.find(st,end))
				hits_minus = len(minus_ranges.find(st,end))
				hits_plus_rpkm = hits_plus*1000000000.0/(size*denominator)
				hits_minus_rpkm = hits_minus*1000000000.0/(size*denominator)
				print >>OUT, '\t'.join(['%s','%d','%d','%s','%d','%s','%d','%d','%.3f','%.3f']) % (chrom,st,end,geneName + "_exon_" + str(exonNum),0,gstrand,hits_plus,hits_minus,hits_plus_rpkm,hits_minus_rpkm)
				if gstrand == '-':exonNum -= 1
				elif gstrand == '+':exonNum += 1
				mRNA_plus_hits += hits_plus
				mRNA_minus_hits += hits_minus
				mRNA_length += size
			mRNA_plus_rpkm = mRNA_plus_hits*1000000000.0/(mRNA_length*denominator)
			mRNA_minus_rpkm = mRNA_minus_hits*1000000000.0/(mRNA_length*denominator)
			print >>OUT, '\t'.join(['%s','%d','%d','%s','%d','%s','%d','%d','%.3f','%.3f']) % (chrom,tx_start,tx_end,geneName + "_mRNA",0,gstrand,mRNA_plus_hits,mRNA_minus_hits,mRNA_plus_rpkm,mRNA_minus_rpkm)
		elif (options.strand_rule is None):	#this is NOT strand specific
			if gstrand == '-':
				intronNum=len(intron_starts)	
				exonNum=len(exon_starts)
			elif gstrand == '+':
				intronNum=1
				exonNum=1
			#assign reads to intron regions
			for st,end in zip(intron_starts,intron_ends):
				if end >st:
					size = end - st
				elif end == st:
					size = 1						
				hits = len(unstrand_ranges.find(st,end))
				hits_rpkm = hits*1000000000.0/(size*denominator)
				print >>OUT, '\t'.join(['%s','%d','%d','%s','%d','%s','%d','%.3f']) % (chrom,st,end,geneName + "_intron_" + str(intronNum),0,gstrand,hits,hits_rpkm)
				if gstrand == '-':intronNum -= 1
				elif gstrand == '+':intronNum +=1
			#assign reads to exon regions
			for st,end in zip(exon_starts,exon_ends):
				if end >st:
					size = end - st
				elif end == st:
					size = 1						
				hits = len(unstrand_ranges.find(st,end))
				hits_rpkm = hits*1000000000.0/(size*denominator)				
				print >>OUT, '\t'.join(['%s','%d','%d','%s','%d','%s','%d','%.3f']) % (chrom,st,end,geneName + "_exon_" + str(exonNum),0,gstrand,hits,hits_rpkm)
				if gstrand == '-':exonNum -= 1
				elif gstrand == '+':exonNum += 1
				mRNA_hits += hits
				mRNA_length += size
			mRNA_rpkm = mRNA_hits*1000000000.0/(mRNA_length*denominator)
			print >>OUT, '\t'.join(['%s','%d','%d','%s','%d','%s','%d','%.3f']) % (chrom,tx_start,tx_end,geneName + "_mRNA",0,gstrand,mRNA_hits,mRNA_rpkm)
		
		gene_finished +=1
		print >>sys.stderr, " %d transcripts finished\r" % (gene_finished),


if __name__ == '__main__':
	main()
