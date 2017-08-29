#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2015.01.27
Usage:	python split_seqs_chimera_check.py seqs.fna output_dir split_part
'''

import sys,os,glob,subprocess

def usage():
    print """
			Usage:	python split_seqs_chimera_check.py seqs.fna output_dir split_part
			output_dir must with a slash
		"""
		
def __main__():

	try:
		seqs = sys.argv[1]
		output_dir = sys.argv[2]
		split_part = sys.argv[3]
	except:
		usage()
		sys.exit(1)
	subprocess.call( 'python ~/scripts/split_fasta_liuchong.py -i  ' + seqs + ' -p ' + split_part + ' -o ' + output_dir, stdout=subprocess.PIPE, shell=True )

	seq_files = os.listdir( output_dir )
	#identify_chimeric_seqs.py -i split_output/seqs.fna -m usearch61 -o usearch_checked_chimeras -r $reference_seqs
	#filter_fasta.py -f split_output/seqs.fna -o chimera_filtered_seqs.fna -s usearch_checked_chimeras/chimeras.txt -n
	for seq_file in seq_files:
		seq_file_ul = output_dir + seq_file
		print "=============Begin to check file %s====================" % seq_file_ul
		chimera_check_output_dir = output_dir + 'usearch_checked_chimeras_' + seq_file
		chimera_ul = chimera_check_output_dir + '/chimeras.txt'
		chimera_filtered_seqs = seq_file_ul + '_filtered.fna'
		subprocess.call( 'identify_chimeric_seqs.py -i ' + seq_file_ul + ' -m usearch61 -o ' + chimera_check_output_dir + ' -r $reference_seqs', stdout=subprocess.PIPE, shell=True )
		subprocess.call( 'filter_fasta.py -f ' + seq_file_ul + ' -o ' + chimera_filtered_seqs + ' -s ' + chimera_ul + ' -n', stdout=subprocess.PIPE, shell=True )
	subprocess.call( 'cat ' + output_dir + '*_filtered.fna > chimera_filtered_seqs.fna', stdout=subprocess.PIPE, shell=True )

if __name__ == "__main__": __main__()	