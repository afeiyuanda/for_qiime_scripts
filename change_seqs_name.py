#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.2.12
Usage:	python change_seqs_name.py split_library_output/seqs.fna mapping.txt new_seqs.fna
'''

import sys

def usage():
    print """
			Usage:	python change_seqs_name.py split_library_output/seqs.fna mapping.txt new_seqs.fna
			-----------------seqs.fna------------------------------------
			>C1.17_1 HWI-M02561:59:000000000-ABTLD:1:1101:18117:2414 orig_bc=ATACACTG new_bc=ATACACTG bc_diffs=0
			TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACCTCTTTCAGTAGGGACGAAGCGCAAGTGACGGTACCTACAGAAGA
			>C1.17_2 HWI-M02561:59:000000000-ABTLD:1:1101:12835:1905 orig_bc=ATACACTG new_bc=ATACACTG bc_diffs=0
			TGAGGAATATTGGACAATGGGCCACAAGCCTGATCCAGCAATTCTGTGTGCACGATGAAGGTCTTCGGATTGTAAAGTGCTTTCAGTTGGGAAGAAGAAA
			
			---------------------mapping.txt------------------------------
			#SampleID	BarcodeSequence	LinkerPrimerSequence	NewSampleID	Description
			C1.17	ATACACTG	CCTACGGGAGGCAGCAG	CKS1	CKS
			C1.21	ATGAGTGT	CCTACGGGAGGCAGCAG	CKS2	CKS
			
			We use this script to change seqs.fna's seqs sample name to mapping file's NewSampleID
		"""

def read_mapping_file( mapping_file_ul, column_flag ):
	#return a dictionary like this: d = {'F':[F1,F2,F3],'M':[M1,M2,M3]}
	f = open( mapping_file_ul )
	lines = f.readlines()
	f.close()
	title_line_list = lines[0].rstrip().split('\t')
	column_flag_index = title_line_list.index( column_flag )
	new_lines = lines[1:] #ingnore the first title line
	group_list = []
	for line in new_lines:
		group = line.rstrip().split('\t')[column_flag_index]
		if group not in group_list:
			group_list.append( group )
	d = {}
	for group in group_list:
		d[group] = []
	print d
	for line in new_lines:
		sampleid = line.split('\t')[0]
		group = line.rstrip().split('\t')[column_flag_index]
		d[group].append( sampleid )
	print d
	return d

def get_seq_sampleid( seq_title ):
	seq_sampleid = seq_title.split()[0].split('_')[0][1:] #>C1.17_1,we only get C1.17
	return seq_sampleid

def write_blocks_into_file( blocks, output_file ):
	outf = open(output_file, 'w')
	for block in blocks:
		outf.writelines( block[0] )
		outf.writelines( block[1] )
	outf.close()
	
def __main__():

	try:
		seqs_file = sys.argv[1]
		mapping_file = sys.argv[2]
		new_seqs_file = sys.argv[3]
	except:
		usage()
		sys.exit(1)
	
	d = read_mapping_file( mapping_file, 'NewSampleID' )
	new_sampleid = d.keys()
	seqs = open(seqs_file).readlines()
	blocks = [seqs[i:i+2] for i in range(0, len(seqs), 2)]
	print "There are %s seqs in this file %s!" % (len(blocks), seqs_file)
	new_blocks = []
	num_flag = 1
	for block in blocks:
		title_line = block[0]
		seq = block[1]
		seq_sampleid = get_seq_sampleid( title_line )
		#print seq_sampleid
		temp_block = []
		for key in d.keys():
			if seq_sampleid in d[key]:
				new_seq_sampleid = key
				new_title_line = '>' + key + '_' + str(num_flag) + ' ' + ' '.join(title_line.split()[1:]) + '\n'
				temp_block.append( new_title_line )
				temp_block.append( seq )
				new_blocks.append(temp_block)
				break
		num_flag = num_flag + 1
	print "We changed %s seqs!" % len(new_blocks)
	print num_flag
	write_blocks_into_file( new_blocks, new_seqs_file )
	
if __name__ == "__main__": __main__()