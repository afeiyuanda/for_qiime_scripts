#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2015.01.15
Usage:	python grep_from_miseq_fastq_20150115.py  input_seqs_ul.fastq  mapping.txt greped.fastq
		edit_output.fastq : 找到从反向引物端开始测序的reads，且以对应的正向引物的反向互补序列结束的reads，然后进行反向互补
		mapping.txt like this:
		#SampleID       BarcodeSequence LinkerPrimerSequence    ReversePrimer   Description
		ND1     TCCCTTGTCTCC    GTGCCAGCMGCCGCGGTAA     GGACTACHVGGGTWTCTAAT    806rcbc0
		ND2     ACGAGACTGATT    GTGCCAGCMGCCGCGGTAA     GGACTACHVGGGTWTCTAAT    806rcbc1
'''
#need cogent/ and fastq.py
import os,sys,re
from fastq import MinimalFastqParser
			
def get_barcodes( mapping_ul ):
	f = open( mapping_ul )
	lines = f.readlines()
	barcodes = []
	for line in lines[1:]:
		temp = line.split( '\t' )
		barcodes.append(temp[1])
	f.close()
	return barcodes

#Reverse complementary sequence
def reverse_com_seq( seq ):
	new_seq = ''.join(["ATCGNFXRYMKSWHBVD"["TAGCNFXRYMKSWHBVD".index(n)] for n in seq[::-1]])
	return new_seq

#reverse a sequence
def reverse_seq( seq ):
	new_seq = seq[::-1]
	return new_seq
	
def usage():
    print """
		Usage:   python grep_from_miseq_fastq_20150115.py  input_seqs_ul.fastq  mapping.txt ouput_dir
		Only have barcodes, no primer seqs.
		#SampleID       BarcodeSequence LinkerPrimerSequence    ReversePrimer   Description
		ND1     TCCCTTGTCTCC    GTGCCAGCMGCCGCGGTAA     GGACTACHVGGGTWTCTAAT    806rcbc0
		ND2     ACGAGACTGATT    GTGCCAGCMGCCGCGGTAA     GGACTACHVGGGTWTCTAAT    806rcbc1
		"""

def __main__():
	try:
		input_seqs_ul = sys.argv[1]
		mapping_ul = sys.argv[2]
		output_dir = sys.argv[3]
	except:
		usage()
		sys.exit(1)
	
	barcodes = get_barcodes( mapping_ul )
	output_fq =  output_dir + '/greped.fq'
	outf = open( output_fq , 'w')
	grep_seqs_log = output_dir + '/grep_seqs_num_log.txt'
	grep_log = open( grep_seqs_log, 'w' )
	
	#initiate an empty dic to record the greped seqs numbers.
	grep_seqs_dic = {}
	for barcode in barcodes:
		grep_seqs_dic[ barcode ] = 0
	
	#no_matched_seqs = 0
	matched_seqs = 0
	total_seqs_num = 0
	
	for seq_id, seq, qual  in MinimalFastqParser(input_seqs_ul, strict=False):
		total_seqs_num = total_seqs_num + 1
		for barcode in barcodes:
			if seq.startswith( barcode ):
				outf.write('@%s\n%s\n+\n%s\n' % (seq_id, seq, qual))
				grep_seqs_dic[ barcode ] = grep_seqs_dic[ barcode ] + 1
				matched_seqs = matched_seqs + 1
				break
			elif seq.endswith( reverse_com_seq(barcode) ):
				new_seq = reverse_com_seq( seq )
				new_qual = reverse_seq( qual )
				outf.write('@%s\n%s\n+\n%s\n' % (seq_id+'|reverse', new_seq, new_qual))
				grep_seqs_dic[ barcode ] = grep_seqs_dic[ barcode ] + 1
				matched_seqs = matched_seqs + 1
				break

	grep_log.write( "There are %s seqs totally!\n" % total_seqs_num )
	grep_log.write( "Matched seqs: %s !\n" % matched_seqs )
	
	lib_name = os.path.basename(input_seqs_ul).split('.')[0]
	for key in grep_seqs_dic:
		grep_log.write( '%s\t%s\t%s\n' % (key, grep_seqs_dic[key], lib_name) )
	grep_log.close()
	outf.close()

if __name__ == "__main__" : __main__()