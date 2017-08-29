#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.15
'''

import sys,os,glob,subprocess,re,argparse

def usage():
    print """
Usage:	python """+sys.argv[0]+ """ 
-indir raw_fq_dir 
-q quality 
-l Discard reads that became shorter than this length 
-outdir /path/to/outdir
		"""
		
def __main__():

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-indir", dest="indir", help="Directory contains all rawdata", required=True);
	newParser.add_argument( "-q", dest="quality", help="quality score", default=20);
	newParser.add_argument( "-l", dest="length", help="Discard reads that became shorter than this length", default=100 );
	newParser.add_argument( "-outdir", dest="outdir", help="Directory contains all clean data", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	indir = argsDict['indir']
	quality = argsDict['quality']
	length = argsDict['length']
	outdir = argsDict['outdir']
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	files1 = glob.glob( indir+'/*.fastq' )
	files2 = glob.glob( indir+'/*.fq' )
	if files1 != []:
		fastq_files = files1
	if files2 != []:
		fastq_files = files2	
	
	for fastq_file in fastq_files:
		cmd = 'trim_galore -q '+str(quality)+' --length '+str(length)+' -o '+outdir+' '+fastq_file
		os.system(cmd)
	
if __name__ == "__main__": __main__()