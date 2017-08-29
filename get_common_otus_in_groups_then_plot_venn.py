#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.30
Usage:	python get_common_otus.py even_otu_table_taxo_name.txt a,b,c
'''

import sys,glob,subprocess,os

def usage():
    print """
			Usage:	python get_common_otus.py even_otu_table_taxo_name.txt mapping.txt group_coloum_name  group_name1,group_name2,group_name3
		"""

def read_mapping_file( mapping_file_ul, column_flag ):
	#return a dictionary like this: d = {'F':[F1,F2,F3],'M':[M1,M2,M3]}
	#Commonly used is Treatment
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

def readfile( file_ul ):
	f = open( file_ul )
	lines = f.readlines()
	new_lines = []
	for line in lines:
		if line.startswith( '#' ):
			continue
		else:
			new_lines.append( line )
	f.close()
	new_lines.insert(0,lines[1][1:]) #row of sample id also has a #
	return new_lines

def judeg_all_zero( line_list, id_list ):
	#if the id_list mapped elements are all 0, return 0, else return 1
	flag = 1
	for id in id_list:
		if line_list[id] == '0.0':
			flag = 0
			break
	return flag		
	
def get_target_group_col_collated( lines, group_dic, group_list ):
	#get_target_group columns,and collated columns from the same group into 1 column, if the columns have one 0, this column marked 0.
	title_line = lines[0]
	samples_list = title_line.split('\t')
	id_dic = {}
	for group in group_list:
		samples = group_dic[group]
		temp_list = []
		for sample in samples:
			temp_list.append( samples_list.index( sample ) )
		id_dic[group] = temp_list
	print "This targeted groups samples' id dic: "
	print id_dic
	new_lines = []
	new_lines.append( '\t'.join(id_dic.keys()) )
	for line in lines[1:]:
		temp_line = []
		line_list = line.strip().split( '\t' )
		#print line_list
		for key in id_dic.keys():
			id_list = id_dic[key]
			flag = judeg_all_zero( line_list, id_list )
			temp_line.append( flag )
		if sum( temp_line )	!= 0:
			new_lines.append( line_list[0] + '\t' + '\t'.join( [str(i) for i in temp_line] ) )
	return new_lines

#line:3951715	0.0	0.0	2.0,the first col is OTU ID
#Then return a new line
def map_otuid_to_samples( line ):
	line_list = line.rstrip().split('\t')
	otuid = line_list[0]
	reads_num_list = line_list[1:]
	new_num_list = []
	for i in reads_num_list:
		if i == '0':
			new_num_list.append('NA')
		else:
			new_num_list.append(otuid)
	return '\t'.join( new_num_list )
	
def __main__():

	try:
		otu_file = sys.argv[1]
		mapping_file = sys.argv[2]
		groups_col_name = sys.argv[3]
		groups = sys.argv[4]
	except:
		usage()
		sys.exit(1)
	group_list = groups.split(',')
	
	otu_lines = readfile( otu_file )
	group_dic = read_mapping_file( mapping_file, groups_col_name )
	
	output_file = 'targeted_groups_formated_otus.txt'
	output = open( output_file, 'w' )
	new_lines = get_target_group_col_collated( otu_lines, group_dic, group_list )
	output.writelines( new_lines[0] + '\n')
	otu_mapping_lines = []
	for line in new_lines[1:]:
		otu_mapping_lines.append( map_otuid_to_samples( line ) )
		
	for line in otu_mapping_lines:
		output.writelines( line + '\n' )
	output.close()
	
	group_num = len( group_list )
	print "There are %s groups!" % group_num
	
	rscript_name = "plot_targeted_group_venn_" + str(group_num) + '.r'
	venn_name = 'venn_diagram_' + str(group_num) + '.tiff'
	Rcmd = open(rscript_name, 'w')
	if group_num == 5:
		Rcmd.writelines(
		"""
library(VennDiagram)
file <- '"""+ output_file +"""'
read.table(file,header=TRUE,sep="\\t") -> venndata
dataf<-as.list(venndata)
new_dataf <- lapply(dataf, function(x) x[!is.na(x)])
mycol = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
venn.plot <- venn.diagram(
new_dataf,
filename = '"""+ venn_name +"""',
col = "black",
fill = mycol,
alpha = 0.50,
cex = c(2.5, 2.5, 2.5, 2.5, 2.5, 2, 0.8, 2, 0.8, 2, 0.8, 2, 0.8,2, 0.8, 2, 0.55, 2, 0.55, 2, 0.55, 2, 0.55, 2, 0.55, 2, 2, 2, 2, 2, 2.5),
cat.col = mycol,
cat.cex = 1,
cat.fontface = "bold",
margin = 0.05)
		""")
	else:
		Rcmd.writelines(
		"""
library(VennDiagram)
file <- '"""+ output_file +"""'
read.table(file,header=TRUE,sep="\\t") -> venndata
dataf<-as.list(venndata)
new_dataf <- lapply(dataf, function(x) x[!is.na(x)])
mycol = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
venn.plot <- venn.diagram(
new_dataf,
filename = '"""+ venn_name +"""',
col = "black",
fill = mycol[1:"""+str(group_num)+"""],
alpha = 0.50,
#cex = c(2.5, 2.5, 2.5, 2.5, 2.5, 2, 0.8, 2, 0.8, 2, 0.8, 2, 0.8,2, 0.8, 2, 0.55, 2, 0.55, 2, 0.55, 2, 0.55, 2, 0.55, 2, 2, 2, 2, 2, 2.5),
cat.col = mycol[1:"""+str(group_num)+"""],
cat.cex = 1,
cat.fontface = "bold",
margin = 0.05)
		""")
		
	Rcmd.close()
	subprocess.call( "R CMD BATCH " + rscript_name, stdout=subprocess.PIPE, shell=True )
	
if __name__ == "__main__": __main__()