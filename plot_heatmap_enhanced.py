#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.02
Usage:	python plot_heatmap.py sorted_otu_table_L6.txt top_num cluster_cols_flag
'''

import sys,subprocess,os

def usage():
    print """
			Usage:	python plot_heatmap.py otu_table_taxo_header_even.txt  top_num  cluster_cols_flag
			Reads in very sample should be normalized!
				otu_table_file = otu_table_taxo_header_even.txt
				top_num = an integer
				cluster_cols_flag = TRUE or FALSE, must be capital!
				row_names otu or taxonomy
		"""

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

def float_list( list ):
	return [float(i) for i in list]
	
def sum_list_convert2str( list ):
	return str( sum(float_list( list )) )

def lines2array( lines_block ):
	#change mulitiple lines into array
	arr = []
	for line in lines_block:
		line_list = line.rstrip().split( '\t' ) # the last line is taxonomy
		otu_seqs_num = sum_list_convert2str( line_list[1:-1] )
		line_list.append( otu_seqs_num )
		arr.append( line_list )
	return arr

#compare several lists' last element, is otu_seqs_num
def cmpl(x,y):
	xtotal = float( x[-1] )
	ytotal = float( y[-1] )
	if xtotal < ytotal:
		return 1
	elif xtotal > ytotal:
		return -1
	else:
		return 0

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
		
def __main__():

	try:
		otu_table_file = sys.argv[1]
		top_num = sys.argv[2]
		cluster_cols_flag = sys.argv[3]
		row_names = sys.argv[4]
	except:
		usage()
		sys.exit(1)
	output_dir = get_dirname(otu_table_file)
	otu_table = get_basename(otu_table_file)
	output = output_dir + '/' + otu_table.split('.')[0] + '_top' + top_num + '.txt'
	lines = readfile( otu_table_file )
	row_names_line = '\t'.join( lines[0].split( '\t' )[:-1] ) + '\n'
	otu_lines = lines[1:]
	lines_arr = lines2array( otu_lines )
	lines_arr.sort( cmp = cmpl )
	outf = open( output, 'w' )
	outf.writelines( row_names_line )
	if row_names == 'otu':
		for lines_list in lines_arr[:int(top_num)]:
			outf.writelines( 'OTU'+'\t'.join(lines_list[:-2]) + '\n' )
		heatmap_name = output_dir+"/heatmap" + '_top' + top_num + '_otu.png'
	else:
		for lines_list in lines_arr[:int(top_num)]:
			outf.writelines( 'OTU'+lines_list[0]+':'+lines_list[-2]+'\t'+'\t'.join(lines_list[1:-2]) + '\n' )
		heatmap_name = output_dir+"/heatmap" + '_top' + top_num + '_taxonomy.png'

	outf.close()
	rscript_name = output_dir+"/plot_heatmap_" + 'top' + top_num + '.r'
	#heatmap_name = output_dir+"/heatmap" + '_top' + top_num + '.png'
	Rcmd = open(rscript_name, 'w')
	Rcmd.writelines(
	"""
library(pheatmap)
otus <- read.table('"""+ output +"""',sep='\\t', header=T, row.names=1, check.names=F )
otu_matrix <- data.matrix(otus)
pheatmap(log(otu_matrix+0.01,2),cluster_cols="""+ cluster_cols_flag +""", cellheight=4, cellwidth =10, color=colorRampPalette(c("#1409a4", "#81ebe5", "#3df134", "#f4f70c","#fb080e"))(50), fontsize=7, fontsize_row=4, border_color = "NA",filename = '"""+ heatmap_name +"""')
	""")
	Rcmd.close()
	subprocess.call( "R CMD BATCH " + rscript_name, stdout=subprocess.PIPE, shell=True )
	#subprocess.call( "convert " + heatmap_name + ' '+output_dir+"/heatmap" + '_top' + top_num + '.png', stdout=subprocess.PIPE, shell=True )
	subprocess.call( "rm -rf *.r *.Rout .RData Rplots.pdf", stdout=subprocess.PIPE, shell=True )
	
if __name__ == "__main__": __main__()


	#lables <- seq(from=ceiling(min(otu_matrix)), to=ceiling(max(otu_matrix)), by=1)
	#breaks<-seq(from=ceiling(min(otu_matrix)), to=ceiling(max(otu_matrix)), by=1)
	#pheatmap(otu_matrix, cellheight = 6,cellwidth = 12,treeheight_col=8,treeheight_row=30,color = colorRampPalette(c("green", "black", "red"))(10), fontsize=7, fontsize_row=4, border_color = "NA",filename = "heatmap.bmp", legend_labels=lables, legend_breaks=breaks)