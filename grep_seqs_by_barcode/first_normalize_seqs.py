#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2015.01.04
'''

import sys,subprocess,re
from fasta import MinimalFastaParser

def usage():
    print """
			Usage:	python first_normalize_seqs.py final_seqs.fna mapping_file normalize_num outfile
		"""	
def read_map_file(map_file):
	sample_names = []
	for line in open(map_file):
		if line.startswith('#') or line.strip()=='':
			pass
		else:
			sample_names.append(line.split('\t')[0])
	return sample_names
		
def __main__():

	try:
		seq_file = sys.argv[1]
		mapping_file = sys.argv[2]
		norm = sys.argv[3]
		outfile = sys.argv[4]
	except:
		usage()
		sys.exit(1)
	
	sample_names = read_map_file(mapping_file)
	sample_seqs_dic = dict(zip(sample_names, [0]*len(sample_names)))
	outf = open(outfile, 'w')
	for seq_id, seq in MinimalFastaParser(open(seq_file), strict=False):
		seq_id = re.split('\s+', seq_id)[0]
		sample_name = seq_id.split('_')[0]
		if sample_seqs_dic[sample_name] < int(norm):
			sample_seqs_dic[sample_name] += 1
			outf.write('>'+sample_name+'_'+str(sample_seqs_dic[sample_name])+'\n'+seq+'\n')
		else:
			pass
	print sample_seqs_dic	
		
if __name__ == "__main__": __main__()