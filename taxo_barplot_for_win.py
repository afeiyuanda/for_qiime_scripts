#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.29
Usage:	python taxo_barplot.py wf_taxa_summary
'''

import sys,glob,subprocess,os

def usage():
    print """
			Usage:	python taxo_barplot.py wf_taxa_summary
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

def sum_every_column( lines_block ):
	if len( lines_block ) == 1:
		return lines_block[0]
	else:
		arr = lines2array( lines_block )
		trans_arr = trans_array( arr )
		sum_list = [ sum(i) for i in trans_arr ]
		sum_str_list = [ str(i) for i in sum_list ]
		sum_line = '\t'.join( sum_str_list )
		final_sum_line = 'Unknown' + '\t' + sum_line + '\n'
	return final_sum_line
	
def change_line_taxa_flag(line):
	#change k__Bacteria;Other to Other,only keep the last word
	unknown_list = ['p__','c__','o__','f__','f__','Other']
	line_list = line.rstrip().split('\t')
	if line_list[0].split(';')[-1] in unknown_list:
		new_line = 'Unknown' + '\t' + '\t'.join(line_list[1:]) + '\n'
	else:
		new_line = line_list[0].split(';')[-1] + '\t' + '\t'.join(line_list[1:]) + '\n'
	return new_line
	
def generate2lines_blocks( lines ):
	#1 block contains unknown taxa lines;
	#2 block contains the other lines
	known_taxa_blocks = []
	unknown_taxa_blocks = []
	unknown_list = ['p__','c__','o__','f__','g__','Other']
	for line in lines:
		if line.rstrip().split('\t')[0].split(';')[-1] in unknown_list:
			unknown_taxa_blocks.append( change_line_taxa_flag(line) )
		else:
			known_taxa_blocks.append( change_line_taxa_flag(line) )
	lines_blocks = []
	lines_blocks.append( known_taxa_blocks )
	lines_blocks.append( unknown_taxa_blocks )
	return lines_blocks

		
def __main__():

	try:
		wf_taxa_summary = sys.argv[1]
	except:
		usage()
		sys.exit(1)
	taxa_file_list = glob.glob( wf_taxa_summary + '/*L*.txt')
	for taxa_file in taxa_file_list:
		level = taxa_file.split('.')[0].split('_')[-1]
		subprocess.call( "mkdir " + level, stdout=subprocess.PIPE, shell=True )
		new_taxa_file_dir = level
		new_taxa_file_ul = level + '/' + level + '.txt'
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
		unknown_taxa_line = sum_every_column( unknown_taxa_block )
		new_taxa_file.writelines( unknown_taxa_line )
		new_taxa_file.close()
		if sample_number <= 10:
			os.chdir(new_taxa_file_dir)
			subprocess.call( "perl E:\\scripts\\R_scripts\\bar_pie_less_10sample.pl -i " + level + '.txt -pie F -bw 4 -bh 5', stdout=subprocess.PIPE, shell=True )
			os.chdir('../')
		else:
			os.chdir(new_taxa_file_dir)
			subprocess.call( "perl E:\\scripts\\R_scripts\\bar_pie.pl -i " + level + '.txt -pie F -bw 4 -bh 5', stdout=subprocess.PIPE, shell=True )
			os.chdir('../')
	
if __name__ == "__main__": __main__()