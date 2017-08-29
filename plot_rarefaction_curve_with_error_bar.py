#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.30
Usage:	python average_collated_alpha_index.py collated_alpha_dir
'''

import sys,glob,subprocess,os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def usage():
    print """
			Usage:	python boxplot_rarefaction_curve.py collated_alpha_dir mapping_file group_flag  legend_pos
			outputs also located in collated_alpha_dir, but with suffix .averaged
			legend pos:¡±bottom¡±,¡±bottomleft¡±,¡±left¡±,¡±topleft¡±,¡±top¡±,¡±topright¡±,¡±right¡±,¡±bottomright¡±,¡±center¡±
		"""

def read_mapping_file( mapping_file_ul, column_flag ):
	#return a dictionary like this: d = {'F':[F1,F2,F3],'M':[M1,M2,M3]}
	#column_flag starts from 0
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

def lines2array( lines_block ):
	#you'd better don't use the same variable name as in the main function
	#change mulitiple lines into array
	arr = []
	for line in lines_block:
		line_list = line.rstrip().split( '\t' )
		arr.append( line_list )
	return arr
	
def trans_array( arr ):
	rows = len( arr )
	columns = len( arr[0] )
	trans_arr = [[r[col]for r in arr] for col in xrange(columns)] 
	return trans_arr

def get_dirname(ul):
	if ul.endswith('/'):
		return os.path.dirname(ul.rstrip('/'))
	else:
		return os.path.dirname(ul)

def get_basename(ul):
	if ul.endswith('/'):
		return os.path.basename(ul.rstrip('/'))
	else:
		return os.path.basename(ul)

def list_all_NA(l):
	flag = True
	for i in l:
		if i != 'NA':
			flag =False
			break
	return flag
		
def __main__():

	try:
		collated_alpha_dir = sys.argv[1]
		mapping_file = sys.argv[2]
		group_flag = sys.argv[3]
		#legend_pos = sys.argv[4]
	except:
		usage()
		sys.exit(1)
	
	group_sample_dic = read_mapping_file( mapping_file, group_flag )
	collated_files = glob.glob( collated_alpha_dir + '/*.txt.averaged' )#the files contain the dir path already.
	for file in collated_files:
		fig = plt.figure()
		f = open( file )
		lines = f.readlines()
		f.close()
		new_group_dic = {}
		for key in group_sample_dic.keys():
			new_group_dic[key] = []
		
		for key in group_sample_dic.keys():
			title_line_list = lines[0].strip().split('\t')
			for sample in group_sample_dic[key]:
				new_group_dic[key].append(title_line_list.index(sample))
		print new_group_dic
		ylab = file.split('.')[0].split('/')[-1]
		output_img= collated_alpha_dir+'/'+ylab + '_errorbar.png'
		output_pdf= collated_alpha_dir+'/'+ylab + '_errorbar.pdf'
		seqs_num_list = []
		for line in lines[1:]:
			seqs_num_list.append(int(line.strip().split('\t')[0]))
		
		x_num = 0
		line_handles = []
		for group in new_group_dic.keys():
			samples_values = []
			for line in lines[1:]:
				tmp_values = []
				for index in new_group_dic[group]:
					tmp_values.append(line.strip().split('\t')[index])
				samples_values.append(tmp_values)
			samples_values_not_empty = []
			for l in samples_values:
				if list_all_NA(l)==False:
					samples_values_not_empty.append(l)
			average_values = []
			standard_errors = []
			for l in samples_values_not_empty:
				tmp_values = []
				for i in l:
					if i != 'NA':
						tmp_values.append(float(i))
				atmp_values = np.array(tmp_values)
				average_values.append(np.mean(atmp_values))
				standard_errors.append(np.std(atmp_values)/np.sqrt(len(atmp_values)-1))
			#print average_values
			#print standard_errors
			if len(average_values)>x_num:
				x_num = len(average_values)
			line_name = 'line_'+group
			line_name, a, b = plt.errorbar(np.array(seqs_num_list[:len(average_values)]), np.array(average_values), yerr=np.array(standard_errors))
			line_handles.append(line_name)
		plt.xlabel('Number of Seqs Sampled')
		plt.ylabel(ylab)
		groups = new_group_dic.keys()
		plt.legend(line_handles,groups, loc=0)
		plt.savefig(output_img)
		#pp = PdfPages(output_pdf)
		#pp.savefig(fig)
		#pp.close()
	
	
if __name__ == "__main__": __main__()