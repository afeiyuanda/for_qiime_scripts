#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.2.12
Usage:	python txt_otu_table2LEfSe_format.py sorted_otu_table_L6.txt output.txt
'''
import sys,os,glob,subprocess,re,argparse
from grep_seqs_by_barcode.fastq import MinimalFastqParser

def usage():
    print """
Usage:	python """+sys.argv[0]+ """ 
-fastq raw_fq 
-mapping /path/to/mapping.txt  
-outdir /path/to/outdir
		"""

#Reverse complementary sequence
def reverse_com_seq( seq ):
	new_seq = ''.join(["ATCGNFXRYMKSWHBVD"["TAGCNFXRYMKSWHBVD".index(n)] for n in seq[::-1]])
	return new_seq
	
def __main__():
	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-fastq", dest="fastq", help="single end fastq format rawdata", required=True);
	newParser.add_argument( "-mapping", dest="mapping", help="mapping file", required=True);
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	fastq = argsDict['fastq']
	mapping = argsDict['mapping']
	outdir = argsDict['outdir']
	subprocess.call( 'dos2unix '+mapping, stdout=subprocess.PIPE, shell=True )
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	primer_index_dic = {}
	sample_primer_dic = {}
	sample_primer_names_list = ['#SampleID','BarcodeSequence','LinkerPrimerSequence','ReversePrimer']
	for line in open(mapping):
		if line.startswith('#'):
			title_line_list = line.strip().split('\t')
			for item in sample_primer_names_list:
				primer_index_dic[item] = title_line_list.index(item)
		else:
			line_list = line.strip().split('\t')
			for item in ['BarcodeSequence','LinkerPrimerSequence','ReversePrimer']:
				sample_primer_dic.setdefault(line_list[primer_index_dic['#SampleID']], []).append(line_list[primer_index_dic[item]])
	print sample_primer_dic
	
	total_seqs_num = 0
	outf = open(outdir + '/greped.fq', 'w')
	log1 = open(outdir + '/grep_log.txt', 'w')
	seq_count_dic = {}
	for sample in sample_primer_dic.keys():
		seq_count_dic[sample] = 0
		barcode_len = len(sample_primer_dic[sample][0])
	try:
		for seq_id, seq, qual  in MinimalFastqParser(fastq, strict=False):
			total_seqs_num = total_seqs_num + 1
			for sample in sample_primer_dic.keys():
				#引物方向均为5 -> 3
				forward_primer = sample_primer_dic[sample][0]
				reverse_primer1 = sample_primer_dic[sample][2]
				reverse_primer2 = sample_primer_dic[sample][0] + sample_primer_dic[sample][2]
				if seq.startswith(forward_primer):
					outf.write('@%s\n%s\n+\n%s\n' % (seq_id, seq, qual))
					seq_count_dic[sample] += 1
					break
				elif seq.startswith(reverse_primer1) and reverse_com_seq(sample_primer_dic[sample][0]) in seq:
					new_seq = sample_primer_dic[sample][0] + reverse_com_seq(seq).split(sample_primer_dic[sample][0])[1]
					new_qual = qual[0:len(new_seq)]
					outf.write('@%s\n%s\n+\n%s\n' % (seq_id, new_seq, new_qual))
					seq_count_dic[sample] += 1
					break
				elif seq.startswith(reverse_primer2):
					try:
						new_seq = forward_primer + reverse_com_seq(seq).split(sample_primer_dic[sample][1])[1]
					except:
						new_seq = forward_primer + reverse_com_seq(seq)
					new_qual = qual[0:len(new_seq)]
					outf.write('@%s\n%s\n+\n%s\n' % (seq_id, new_seq, new_qual))
					seq_count_dic[sample] += 1
					break
	except:
		pass
	outf.close()
	print seq_count_dic
	for sample in seq_count_dic.keys():
		log1.write(sample + '\t' + str(seq_count_dic[sample]) + '\n')
	log1.close()
	subprocess.call( 'perl '+sys.path[0]+'/FastqToFasta.pl -i ' + outdir+'/greped.fq', stdout=subprocess.PIPE, shell=True )
	
	subprocess.call( 'split_libraries.py -l 100 -f '+outdir+'/greped.fq_fasta'+' -b '+str(barcode_len)+' -m ' + mapping + ' -o ' + outdir, stdout=subprocess.PIPE, shell=True )
	
	subprocess.call( 'identify_chimeric_seqs.py -i '+outdir+'/seqs.fna -m usearch61 -o '+outdir+'/usearch_checked_chimeras -r $reference_seqs', stdout=subprocess.PIPE, shell=True )
	#filter_fasta.py -f split_output/seqs.fna -o chimera_filtered_seqs.fna -s usearch_checked_chimeras/chimeras.txt -n
	subprocess.call( 'filter_fasta.py -f '+outdir+'/seqs.fna -o '+outdir+'/usearch_checked_chimeras/chimera_filtered_seqs.fna -s '+\
		outdir+'/usearch_checked_chimeras/chimeras.txt -n', stdout=subprocess.PIPE, shell=True )
	os.system('mv '+outdir+'/usearch_checked_chimeras/chimera_filtered_seqs.fna ' + outdir+'/final_seqs.fna' )
	
if __name__ == "__main__": __main__()