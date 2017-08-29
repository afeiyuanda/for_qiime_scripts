#!/usr/bin/python  
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2013.1.15, modified 2014.05.19
Usage:	Usage:	python biopython_xml.py xml_fomat.blast_out  one_line.blast_out
		Firstly, change xml blast output into flat blast output;
		Seconde, delete the duplication rows by the first column.

Linux shell delete the duplication rows by the first column.
awk -F '\t' '$1 == cust_id{} $1 != cust_id{cust_id = $1;print $0}' test.txt  > end_result.txt

alignment.accession  alignment.hit_id     alignment.length
alignment.hit_def    alignment.hsps       alignment.title

hsp.align_length    hsp.frame           hsp.match           hsp.query           hsp.sbjct           hsp.score           
hsp.bits            hsp.gaps            hsp.num_alignments  hsp.query_end       hsp.sbjct_end       hsp.strand          
hsp.expect          hsp.identities      hsp.positives       hsp.query_start     hsp.sbjct_start

blast_record.alignments
blast_record.application
blast_record.blast_cutoff
blast_record.database
blast_record.database_length
blast_record.database_letters
blast_record.database_name
blast_record.database_sequences
blast_record.date
blast_record.descriptions
blast_record.dropoff_1st_pass
blast_record.effective_database_length
blast_record.effective_hsp_length
blast_record.effective_query_length
blast_record.effective_search_space
blast_record.effective_search_space_used
blast_record.expect
blast_record.filter
blast_record.frameshift
blast_record.gap_penalties
blast_record.gap_trigger
blast_record.gap_x_dropoff
blast_record.gap_x_dropoff_final
blast_record.gapped
blast_record.hsps_gapped
blast_record.hsps_no_gap
blast_record.hsps_prelim_gapped
blast_record.hsps_prelim_gapped_attemped
blast_record.ka_params
blast_record.ka_params_gap
blast_record.matrix
blast_record.multiple_alignment
blast_record.num_good_extends
blast_record.num_hits
blast_record.num_letters_in_database
blast_record.num_seqs_better_e
blast_record.num_sequences
blast_record.num_sequences_in_database
blast_record.posted_date
blast_record.query
blast_record.query_id
blast_record.query_length
blast_record.query_letters
blast_record.reference
blast_record.sc_match
blast_record.sc_mismatch
blast_record.threshold
blast_record.version
blast_record.window_size

'''
from Bio.Blast import NCBIXML
import os,sys,subprocess,re

def usage():
    print """
			Usage:	python biopython_xml.py xml_fomat.blast_out
			Firstly, change xml blast output into flat blast output;
			Seconde, delete the duplication rows by the first column;
			Final, output file's  suffix is .flat_deduplicate
		"""

def get_taxa_info(desc):
	p = re.compile(r'\[(.*?)\]')
	m = p.findall(desc)
	length_tax = ''
	for i in m:
		if len(i)>len(length_tax):
			length_tax = i
	if length_tax.startswith('['):
		length_tax = length_tax[1:]
	return length_tax
		
def __main__():
	
	try:
		blast_out = sys.argv[1]
		cutoff = sys.argv[2]
		#output_file = sys.argv[2]
	except:
		usage()
		sys.exit(1)
	
	output_file = blast_out.split('.')[0] + '.flat'
	final_output = output_file + '_deduplicate'	
	if os.path.exists(final_output):
		pass
	else:
		result_handle = open( blast_out )
		blast_records = NCBIXML.parse( result_handle )
		#blast_record = blast_records.next()
		output = open( output_file, 'w' )
		output.writelines( 'Query\tBest Hit\tQ.start\tQ.end\tS.start\tS.end\tAlignment Length\tIdentity\tBit Score\tE-value\tTaxonomy\n' )
		try:
			for blast_record in blast_records:
				for alignment in blast_record.alignments:
					#print alignment.hsps[0] 输出序列比对信息
					#print alignment 输出alignment.hit_def 和length
					for hsp in alignment.hsps:
						output.writelines( blast_record.query + '\t' \
						+ alignment.hit_id + '\t' \
						+ str(hsp.query_start) + '\t' \
						+ str(hsp.query_end) + '\t' \
						+ str(hsp.sbjct_start) + '\t' \
						+ str(hsp.sbjct_end) + '\t' \
						+ str(hsp.align_length) + '\t' \
						#+ str(hsp.identities) + '\t' \
						+ str( round(float(hsp.identities)/hsp.align_length, 4) ) +'\t'\
						+ str(hsp.bits) + '\t' \
						+ str(hsp.expect) + '\t' + get_taxa_info(alignment.hit_def) + '\n' )
						#+ str(hsp.expect) + '\t' + alignment.hit_def.split('|')[4].split('[')[0].strip() + '\n' )
		except:
			print 'error'
	
		print "Delete the duplication rows by the first column, and write the output into %s!" % final_output
		subprocess.call( 'python '+sys.path[0]+'/blast_fomat8_file_dereplication.py '+output_file, stdout=subprocess.PIPE, shell=True )	
	work_dir = os.path.dirname(blast_out)
	taxa_output = open(work_dir+'/taxa_assign.txt', 'w')
	for line in open(final_output):
		if line.startswith('Query'):
			pass
		else:
			line_list = line.strip().split('\t')
			if float(line_list[7]) > float(cutoff):
				try:
					taxa_output.write(re.split('\s+',line_list[0])[0]+'\t'+line_list[10]+'\t'+line_list[7]+'\n')
				except:
					print 'This taxa line do not have taxa name!'
			else:
				taxa_output.write(re.split('\s+',line_list[0])[0]+'\tUnknown\t'+line_list[7]+'\n')
	taxa_output.close()

if __name__ == "__main__" : __main__()
