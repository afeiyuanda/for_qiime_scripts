#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.15
'''

import sys,os,glob,subprocess,re,argparse

def usage():
    print """
Usage:	python """+sys.argv[0]
	
def __main__():

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-biom", dest="biom", help="otu_table_mc2_w_tax_no_pynast_failures.biom", required=True);
	newParser.add_argument( "-mapping", dest="mapping", help="mapping file");
	newParser.add_argument( "-filter", dest="filter", help="filtered samples");
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	biom = argsDict['biom']
	mapping = argsDict['mapping'] 
	filter = argsDict['filter']
	outdir = argsDict['outdir']
	
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	samples_list = []
	for line in open(mapping):
		if line.startswith('#'):
			pass
		else:
			samples_list.append(line.split('\t')[0])
	filtered_sample = [i.strip() for i in filter.strip().split(',')]
	kept_samples = []
	for s in samples_list:
		if s not in filtered_sample:
			kept_samples.append(s)
	kept = open(outdir+'/kept_samples.txt', 'w')
	kept.write('\n'.join(kept_samples))
	kept.close()
	#filter_samples_from_otu_table.py -i ../picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom -m ../mapping.txt --output_mapping_fp mapping_filter.txt  --sample_id_fp sample_kept.txt -o filter.biom
	filter_cmd = 'filter_samples_from_otu_table.py -i '+biom+' -m '+mapping+' --output_mapping_fp '+outdir+'/mapping_filtered.txt  --sample_id_fp '+outdir+'/kept_samples.txt -o '+outdir+'/filtered.biom'
	print filter_cmd
	subprocess.call( filter_cmd, stdout=subprocess.PIPE, shell=True )
		
if __name__ == "__main__": __main__()