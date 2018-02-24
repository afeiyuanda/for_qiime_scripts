#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2015.01.30
Modify:	2016.07.28,当组只包含一个样品时，进行了剔除，然后再进行组与组间差异分析
Usage:	python taxonomy_anova.py taxa_summary_even/ mapping_file group_flag
'''

import sys,copy,subprocess,glob

def usage():
    print """
			Usage:	python taxonomy_anova.py taxa_summary_even/ mapping_file group_flag
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
	
	for key in groups_dic.keys():
		if len(groups_dic[key]) == 1:
			del groups_dic[key]
	
	min_samples = 10000 #keep the min sample NO. in all groups
	max_samples = 0
	for key in groups_dic.keys():
		if len( groups_dic[key] ) < min_samples:
			min_samples = len( groups_dic[key] )
		if len( groups_dic[key] ) > max_samples:
			max_samples = len( groups_dic[key] )
			
	result_list = []
	result_list.append( groups_dic )
	#print groups_dic
	result_list.append( empty_groups_dic )
	result_list.append( min_samples )
	result_list.append( max_samples )
	return result_list

def trans_list2array(list2array):
    rows = len(list2array)
    columns = len(list2array[0])
    t_list2array=[[r[col]for r in list2array] for col in xrange(columns)]
    return t_list2array
	
'''
Final temp_taxo_file like this:
        Others  Group
Feed1   0.0133571428571394      Feed
Feed2   0.0134999999999999      Feed
GSF1    0.0130714285714282      GSF
GSF2    0.0118571428571202      GSF
'''	
def generate1taxa_file( title_line, taxa_line , group_list, temp_taxo_file):
	lines_list = []
	lines_list.append( title_line.rstrip().split('\t') )
	lines_list.append( taxa_line.rstrip().split('\t') )
	lines_list.append( group_list )
	trans_arrary = trans_list2array( lines_list )
	f = open( temp_taxo_file, 'w' )
	trans_lines = [ '\t'.join(i) for i in trans_arrary ]
	for line in trans_lines:
		f.writelines( line + '\n' )
	f.close()

def generate_group_by_group_list( group_list ):
	group_by_group_list = []
	for i in range( len(group_list)-1 ):
		j = i + 1
		while (j < len(group_list)):
			group_by_group_list.append( group_list[i] + '-' + group_list[j] )
			j = j + 1
	return group_by_group_list

def get_pvalue_list( group_by_group_list, anova_test_output ):
	anova_lines = open( anova_test_output ).readlines()
	pvalue_list = []
	for group_by_group in group_by_group_list:
		rv_group_by_group = group_by_group.split('-')[1] + '-' + group_by_group.split('-')[0]
		for line in anova_lines:
			if group_by_group in line or rv_group_by_group in line:
				pvalue_list.append( line.rstrip().split(' ')[-1] )
	return pvalue_list	
	
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
	
	groups = groups_dic.keys()
	group_by_group_list = generate_group_by_group_list( groups )
	
	Rcmd = open("anova_test.r", 'w')
	Rcmd.writelines(
	"""
mydata <- read.table('temp_taxo.txt',check.names=FALSE,header=T,row.names=1,sep='\t')
aov_test <- aov(mydata[,1]~as.factor(mydata[,2]))
sink("anova_test_output.txt")
TukeyHSD(aov_test)
sink()
	"""
	)
	Rcmd.close()
	
	level_files = glob.glob( taxa_summary_even + '/L*')
	for level_file in level_files:
		empty_groups_dic = get_groups_info(mapping_file, group_flag)[1]
		#{'M': [], 'D': [], 'F': []}
		original_taxo_file = level_file + '/average.xls'
		print original_taxo_file
		f = open( original_taxo_file )
		lines = f.readlines()
		f.close()
		title_line = lines[0]
		samples = title_line.strip().split('\t')
		new_samples = []
		group_list = ['Group']
		for sample in samples:
			for key in groups_dic.keys():
				#del containing only 1 sample group
				if sample in groups_dic[key] and len(groups_dic[key])>1:
					group_list.append( key )
					new_samples.append(sample)
		new_title_line = '\t'.join(new_samples)
					
		final_test_output = 'anova_test.xls'
		anova_test_result = level_file + '/' + final_test_output
		outputf = open( anova_test_result, 'w' )
		new_title_list = title_line.rstrip().split('\t') + group_by_group_list
		outputf.writelines( '\t'.join(new_title_list) + '\n' )
		
		for taxa_line in lines[1:]:
			generate1taxa_file( new_title_line, taxa_line , group_list, 'temp_taxo.txt')
			subprocess.call( "R CMD BATCH anova_test.r", stdout=subprocess.PIPE, shell=True )
			pvalue_list = get_pvalue_list( group_by_group_list, "anova_test_output.txt" )
			taxa_line_list = taxa_line.rstrip().split('\t')
			new_taxa_line = '\t'.join( taxa_line_list + pvalue_list )
			outputf.writelines( new_taxa_line + '\n' )
		
		outputf.close()
if __name__ == "__main__": __main__()
