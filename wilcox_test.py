#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2015.01.04
Usage:	python wicox_test.py mapping_file/  taxa_summary_even
'''

import sys,copy,subprocess,glob

def usage():
    print """
			Usage:	python wicox_test.py taxa_summary_even/ mapping_file group_flag
			Notice the / !
		"""
		
#This fun return a list,list[0] is groups and its samples formed dic, list[1] is groups empty dic, list[2] is the min sample NO. in all groups
def get_groups_info(mapping_file_path, group_flag):
	f = open( mapping_file_path )
	lines = f.readlines()
	groups_dic = {}
	groups = []
	for line in lines:
		if line.startswith('#'):
			group_flag_index = line.rstrip().split('\t').index( group_flag )
			continue
		else:
			group = line.rstrip().split('\t')[group_flag_index]
			if group not in groups:
				groups.append( group )
	print groups	
	for group in groups:
		groups_dic[group] = []
	empty_groups_dic = copy.deepcopy(groups_dic) #we use copy.deepcopy, intend to make a new dic that will not change when groups_dic change
	
	for group in groups:
		for line in lines:
			if line.rstrip().split('\t')[group_flag_index]==group:
				groups_dic[group].append( line.rstrip().split('\t')[0] )
	min_samples = 10000 #keep the min sample NO. in all groups
	max_samples = 0
	for group in groups:
		if len( groups_dic[group] ) < min_samples:
			min_samples = len( groups_dic[group] )
		if len( groups_dic[group] ) > max_samples:
			max_samples = len( groups_dic[group] )
			
	result_list = []
	result_list.append( groups_dic )
	#print groups_dic
	result_list.append( empty_groups_dic )
	result_list.append( min_samples )
	result_list.append( max_samples )
	return result_list
		
def __main__():

	try:
		taxa_summary_even = sys.argv[1]
		mapping_file = sys.argv[2]
		group_flag = sys.argv[3]
	except:
		usage()
		sys.exit(1)
	groups_dic = get_groups_info(mapping_file, group_flag)[0]
	groups_no = len(groups_dic.keys())
	#{'M': ['M1', 'M2', 'M3'], 'D': ['D1', 'D2', 'D3'], 'F': ['F1', 'F2', 'F3']}
	
	level_files = glob.glob( taxa_summary_even + '/L*')
	for level_file in level_files:
		empty_groups_dic = get_groups_info(mapping_file, group_flag)[1]
		#{'M': [], 'D': [], 'F': []}
		original_taxo_file = level_file + '/average.xls'
		print original_taxo_file
		f = open( original_taxo_file )
		lines = f.readlines()
		samples = lines[0].rstrip().split('\t')
		for group in groups_dic.keys():
			for sample in groups_dic[group]:
				empty_groups_dic[group].append( samples.index( sample ) + 1)
		print empty_groups_dic
		#empty_groups_dic now changed as following:
		#{'M': [7, 8, 9], 'D': [1, 2, 3], 'F': [4, 5, 6]} all add 1
		#sort dic by value
		sample_id_list = sorted(empty_groups_dic.items(), lambda x, y: cmp(x[1], y[1]))
		group_data_list = []
		for group in sample_id_list:
			data_line = 'data' + group[0] + '<-mydata[i,][' + str(min(group[1])) + ':' + str(max(group[1])) + ']\n' 
			group_data_list.append( data_line )
		p_value_list = []
		empty_p_value_list = []
		data_list = []
		for d in sample_id_list:
			data_list.append( d[0] )
		for data in data_list:
			data_index = data_list.index( data )
			i = data_index + 1
			while i < len( data_list ):
				p_value_list.append( 'mydata$Group_' + data + '_' + data_list[i] + '_p_value[i]<-wilcox.test(data.matrix(data' + data + '),data.matrix(data' + data_list[i] + '),alternative="two.sided",exact=FALSE,correct=FALSE)[[3]]\n' )
				empty_p_value_list.append( 'mydata$Group_' + data + '_' + data_list[i] + '_p_value<-c()\n' )
				i = i + 1
			
		Rcmd = open("wilcox_test.r", 'w')
		Rcmd.writelines( 'mydata<-read.table(file="' + original_taxo_file + '",sep="\\t",header=TRUE,check.names=FALSE)\n' )
		for line in empty_p_value_list:
			Rcmd.writelines( line )
		Rcmd.writelines( 'for (i in 1:length(rownames(mydata))){\n' )
		for line in group_data_list:
			Rcmd.writelines( line )
		for line in p_value_list:
			Rcmd.writelines( line )
		Rcmd.writelines('}\n')
		
		final_test_output = 'average_wilcox_test_' + '_'.join( groups_dic.keys() ) + '.xls'
		wilcox_test_result = level_file + '/' + final_test_output
		Rcmd.writelines( 'write.table(mydata,file="' + wilcox_test_result + '",sep="\\t",col.names=TRUE,row.names=FALSE,quote=FALSE)' )
		Rcmd.close()
		subprocess.call( "R CMD BATCH wilcox_test.r", stdout=subprocess.PIPE, shell=True )
		
if __name__ == "__main__": __main__()