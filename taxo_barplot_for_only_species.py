#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.29
Usage:	python taxo_barplot.py wf_taxa_summary output_dir
'''

import sys,glob,subprocess,os,re

def usage():
    print """
			Usage:	python taxo_barplot.py wf_taxa_summary/ output_dir
		"""

def lines2array( lines_block ):
	#you'd better don't use the same variable name as in the main function
	#change mulitiple lines into array
	arr = []
	for line in lines_block:
		line_list = line.rstrip().split( '\t' )[1:]
		line_int_list = []
		for i in line_list:
			line_int_list.append( float(i) )
		arr.append( line_int_list )
	return arr
	
def trans_array( arr ):
	rows = len( arr )
	columns = len( arr[0] )
	trans_arr = [[r[col]for r in arr] for col in xrange(columns)] 
	return trans_arr

def sum_same_taxo_column( lines_block ):
	arr = lines2array( lines_block )
	trans_arr = trans_array( arr )
	sum_list = [ sum(i) for i in trans_arr ]
	sum_str_list = [ str(i) for i in sum_list ]
	sum_line = '\t'.join( sum_str_list )
	final_sum_line = lines_block[0].split('\t')[0] + '\t' + sum_line + '\n'
	return final_sum_line
	
def change_line_taxa_flag(line):
	#change k__Bacteria;Other to Other,only keep the last word
	unknown_list = ['p__','c__','o__','f__','g__','Other','p__unidentified','c__unidentified','o__unidentified','f__unidentified','g__unidentified','p__Incertae_sedis','c__Incertae_sedis','o__Incertae_sedis','f__Incertae_sedis','g__Incertae_sedis','Unknown','None']
	line_list = line.rstrip().split('\t')
	#p = re.compile(r'(\w__)\[(\S+)\]')
	if line_list[0] in unknown_list:
		new_line = 'Unknown' + '\t' + '\t'.join(line_list[1:]) + '\n'
	else:
		new_line = line
	return new_line
	
def generate2lines_blocks( lines ):
	#1 block contains unknown taxa lines;
	#2 block contains the other lines
	known_taxa_blocks = []
	unknown_taxa_blocks = []
	unknown_list = ['p__','c__','o__','f__','g__','Other','p__unidentified','c__unidentified','o__unidentified','f__unidentified','g__unidentified','p__Incertae_sedis','c__Incertae_sedis','o__Incertae_sedis','f__Incertae_sedis','g__Incertae_sedis','Unknown','None']
	known_taxa_dic = {}
	for line in lines:
		if line.rstrip().split('\t')[0] in unknown_list:
			unknown_taxa_blocks.append( change_line_taxa_flag(line) )
		else:
			new_line = change_line_taxa_flag(line)
			taxa_name = new_line.split('\t')[0]
			known_taxa_dic.setdefault(taxa_name,[]).append(new_line)
			
	for taxa_name in known_taxa_dic:
		if len(known_taxa_dic[taxa_name]) > 1:
			known_taxa_blocks.append( sum_same_taxo_column(known_taxa_dic[taxa_name]) )
		else:
			known_taxa_blocks.append( known_taxa_dic[taxa_name][0] )
	lines_blocks = []
	lines_blocks.append( known_taxa_blocks )
	lines_blocks.append( unknown_taxa_blocks )
	return lines_blocks

		
def __main__():

	try:
		wf_taxa_summary = sys.argv[1]
		output_dir = sys.argv[2]
	except:
		usage()
		sys.exit(1)
	taxa_file_list = glob.glob( wf_taxa_summary + '/*L*.txt')
	if os.path.exists( output_dir ):
		print "%s is existing already!" % output_dir
	else:
		os.mkdir( output_dir )
	for taxa_file in taxa_file_list:
		level = taxa_file.split('.')[0].split('_')[-1]
		new_taxa_file_dir = output_dir + '/' + level
		subprocess.call( 'mkdir -p ' + new_taxa_file_dir, stdout=subprocess.PIPE, shell=True )
		new_taxa_file_ul = new_taxa_file_dir + '/' + level + '.txt'
		new_taxa_file = open( new_taxa_file_ul, 'w' )
		
		taxaf = open( taxa_file )
		lines = taxaf.readlines()
		taxaf.close()
		sample_number = len(lines[0].rstrip().split( '\t' ))-1
		
		new_taxa_file.writelines( lines[0] )
		two_blocks = generate2lines_blocks( lines[1:] )
		known_taxa_block = two_blocks[0]
		for l in known_taxa_block:
			new_taxa_file.writelines( l )
			
		unknown_taxa_block = two_blocks[1]
		print unknown_taxa_block
		unknown_taxa_line = sum_same_taxo_column( unknown_taxa_block )
		print unknown_taxa_line
		new_taxa_file.writelines( unknown_taxa_line )
		new_taxa_file.close()
		if sample_number <= 10:
			os.chdir(new_taxa_file_dir)
			subprocess.call( 'perl '+sys.path[0]+'/Rscripts/bar_pie_less_10sample_for_species.pl -i ' + level + '.txt -pie F -bw 4 -bh 4', stdout=subprocess.PIPE, shell=True )
			#subprocess.call( 'rm -rf *.r *.Rout', stdout=subprocess.PIPE, shell=True )
			#os.chdir('../../')
		else:
			os.chdir(new_taxa_file_dir)
			subprocess.call( 'perl '+sys.path[0]+'/Rscripts/bar_pie_for_species.pl -i ' + level + '.txt -pie F -bw 4 -bh 4', stdout=subprocess.PIPE, shell=True )
			#subprocess.call( 'rm -rf *.r *.Rout', stdout=subprocess.PIPE, shell=True )
			#os.chdir('../../')
	
if __name__ == "__main__": __main__()