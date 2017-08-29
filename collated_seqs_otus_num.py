#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.2.12
Usage:	python collated_seqs_otus_num.py mapping_file.txt  Treatment  seqs_per_sample.txt otus_per_sample.txt even_seqs_per_sample.txt even_otus_per_sample.txt total_summary.txt
'''

import sys

def usage():
    print """
			Usage:	python collated_seqs_otus_num.py mapping_file.txt  group_col_name  seqs.txt otus.txt even_seqs.txt even_otus.txt total_summary.txt
			group_col_name: group column name, group_col: starts from 0, identify group information
			seqs_per_sample_file:start_line_flag is 16
			otus_per_sample_file:start_line_flag is 14
		"""

def read_seqs_file(seqs_per_sample_file,start_line_flag):
	#seqs_per_sample_file:start_line_flag is 16
	#otus_per_sample_file:start_line_flag is 14
	f = open( seqs_per_sample_file )
	lines = f.readlines()
	f.close()
	d = {}
	for line in lines[int(start_line_flag):]:
		line_list = line.rstrip().split(' ')
		sampleid = line_list[1][:-1]
		seqs_num = line_list[2]
		d[sampleid] = seqs_num
	print d
	return d
	
def __main__():

	try:
		mapping_file = sys.argv[1]
		group_col_name = sys.argv[2]
		seqs = sys.argv[3]
		otus = sys.argv[4]
		even_seqs = sys.argv[5]
		even_otus = sys.argv[6]
		total_summary = sys.argv[7]
	except:
		usage()
		sys.exit(1)
	
	f = open( mapping_file )
	lines = f.readlines()
	f.close()
	seqs_per_sample_dic = read_seqs_file(seqs,'16')
	otus_per_sample_dic = read_seqs_file(otus,'14')
	even_seqs_per_sample_dic = read_seqs_file(even_seqs,'16')
	even_otus_per_sample_dic = read_seqs_file(even_otus,'14')
	
	otus_num = open(seqs).readlines()[1].rstrip().split(' ')[-1]
	seqs_num = open(seqs).readlines()[2].rstrip().split(' ')[-1]
	even_otus_num = open(even_seqs).readlines()[1].rstrip().split(' ')[-1]
	even_seqs_num = open(even_seqs).readlines()[2].rstrip().split(' ')[-1]
	
	f = open(total_summary, 'w')
	title_line = '\t'.join(['SampleID','Group','SeqsNum','OTUsNum','EvenSeqsNum','EvenOTUsNum'])
	f.writelines(title_line + '\n')
	group_col = lines[0].split('\t').index(group_col_name)
	for line in lines[1:]:
		l = line.rstrip().split('\t')
		sampleid = l[0]
		group = l[int(group_col)]
		try:
			new_line = '\t'.join([sampleid, group, seqs_per_sample_dic[sampleid], otus_per_sample_dic[sampleid], even_seqs_per_sample_dic[sampleid], even_otus_per_sample_dic[sampleid]])
		except:
			new_line = '\t'.join([sampleid, group, seqs_per_sample_dic[sampleid], otus_per_sample_dic[sampleid], '-', '-'])
		f.writelines(new_line + '\n')
	f.writelines( '%s\t\t%s\t%s\t%s\t%s\n' % ('Total', seqs_num, otus_num, even_seqs_num, even_otus_num ) )
	f.close()	
	
if __name__ == "__main__": __main__()