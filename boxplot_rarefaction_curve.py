#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.30
Usage:	python average_collated_alpha_index.py collated_alpha_dir
'''

import sys,glob,subprocess,os

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
	
def __main__():

	try:
		collated_alpha_dir = sys.argv[1]
		mapping_file = sys.argv[2]
		group_flag = sys.argv[3]
		legend_pos = sys.argv[4]
	except:
		usage()
		sys.exit(1)
	
	collated_files = glob.glob( collated_alpha_dir + '/*.txt.averaged' )#the files contain the dir path already.
	for file in collated_files:
		f = open( file )
		lines = f.readlines()
		f.close()
		out_lines = []
		output_ul = file + '.trans'
		outf = open( output_ul, 'w' )
		lines_arr = lines2array( lines )
		t_lines_arr = trans_array( lines_arr )
		for line_arr in t_lines_arr:
			out_line = '\t'.join(line_arr) + '\n'
			out_lines.append( out_line )
			outf.writelines( out_line )
		outf.close()

	group_sample_dic = read_mapping_file( mapping_file, group_flag )
	new_group_dic = {}
	for key in group_sample_dic.keys():
		new_group_dic[key] = []
		
	for key in group_sample_dic.keys():
		for line in out_lines:
			for sample in group_sample_dic[key]:
				if line.startswith( sample ):
					new_group_dic[key].append(out_lines.index(line))
	print new_group_dic
		
	averaged_files = glob.glob( collated_alpha_dir+'/*.averaged' )
	for file in averaged_files:
		Rcmd = open(collated_alpha_dir+"/boxplot_alpha_index.r", 'w')
		ylab = file.split('.')[0]
		output_img= ylab + '.box_plot.png'
		real_ylab = ylab.split('/')[-1]
		boxplot_lines = []
		i = 1
		for group in new_group_dic.keys():
			if i == 1:
				plot_line = """boxplot(data[""" + str(new_group_dic[group][0]) + """:""" + str(new_group_dic[group][-1]) + """,],col=mycol[""" + str(i) + """],xlab="Number of Seqs Sampled",ylab='""" + real_ylab + """',ylim=c(min(data[!is.na(data)]),max(new_data)*(1+0.02)) )\n""" 
				boxplot_lines.append( plot_line )
			else:
				plot_line = """boxplot(data[""" + str(new_group_dic[group][0]) + """:""" + str(new_group_dic[group][-1]) + """,],col=mycol[""" + str(i) + """],xlab="Number of Seqs Sampled",ylab='""" + real_ylab + """' ,ylim=c(min(data[!is.na(data)]),max(new_data)*(1+0.02)) , add=T)\n""" 
				boxplot_lines.append( plot_line )
			i = i + 1
		
		Rcmd.writelines(
		"""
library(RColorBrewer)
mycol <- c( brewer.pal(9,"Set1"),brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(8,"Accent"))
mycol <- rep(mycol, 3)
d<-read.table(file='"""+ file +"""', header=T, row.names=1)
data <- t(d)
new_data <- read.table(file='"""+ file +"""', header=T, row.names=1)
new_data[is.na(new_data)] <- 0
png('"""+ output_img + """',res = 300, width = 3500, height = 3000)
#plot(x=colnames(data), y=x, ylim=c(min(new_data), max(new_data)), xlab="Number of Seqs Sampled",ylab='"""+ real_ylab +"""', type="n", font.lab=1)
""")
	
		for plot_line in boxplot_lines:
			Rcmd.writelines( plot_line )
		
		Rcmd.writelines(
"""legend('""" + legend_pos + """',legend=c('""" + "','".join(new_group_dic.keys()) + """'), col=mycol[1:""" + str(len(new_group_dic)) +"""],lty=1,lwd=2.5,cex=1,ncol=1,horiz=F)
dev.off()
		""")
		Rcmd.close()
		subprocess.call( "R CMD BATCH "+collated_alpha_dir+"/boxplot_alpha_index.r", stdout=subprocess.PIPE, shell=True )
		subprocess.call( "rm -rf *.r *.Rout .RData Rplots.pdf", stdout=subprocess.PIPE, shell=True )
		#subprocess.call( "convert " + output_img + ' ' + ylab + '.png' , stdout=subprocess.PIPE, shell=True )
	
if __name__ == "__main__": __main__()