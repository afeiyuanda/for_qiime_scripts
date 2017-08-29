#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.30
Usage:	python average_collated_alpha_index.py collated_alpha_dir lengend_pos
'''

import sys,glob,subprocess,os

def usage():
    print """
			Usage:	python plot_pca_with_pc_vectors.py pcoa_weighted_unifrac_otu_table_even1000.txt mapping.txt Treatment legend_pos
			Through Treatment column group samples,you can also choose the other col.
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
	
def __main__():

	try:
		pcoa_file = sys.argv[1]
		mapping_file = sys.argv[2]
		group_flag = sys.argv[3]
		lengend_pos = sys.argv[4]
	except:
		usage()
		sys.exit(1)
	
	d = read_mapping_file(mapping_file, group_flag)
	f = open( pcoa_file )
	lines = f.readlines()
	f.close()
	variation_line = lines[-1]
	pc1_var = lines[-1].split('\t')[1]
	pc2_var = lines[-1].split('\t')[2]
	title_line = 'group\t' + lines[0]
	output_name = pcoa_file + '.formated'
	outf = open( output_name, 'w')
	outf.writelines( title_line )
	for line in lines[1:]:
		if line == '\n':
			break
		else:
			sampleid = line.split('\t')[0]
			for key in d.keys():
				if sampleid in d[key]:
					group = key
					break
			new_line = group + '\t' + line
			outf.writelines( new_line )
	outf.close()
	output_img = pcoa_file.split('.')[0] + '.tiff'
	Rcmd = open("plot_pca.r", 'w')	
	Rcmd.writelines(
		"""
library(RColorBrewer)
col_list <- c( brewer.pal(9,"Set1"),brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(8,"Accent"))
sample.bg <- col_list[1:31]

d <- read.table('"""+ output_name +"""',header=T,sep='\\t')
pc1 <- round( """+ pc1_var +""", 2 )
pc2 <- round( """+ pc2_var +""", 2 )
lable <- '% variation explained'
xlab <- paste( 'PC1 (', pc1, lable,')', sep='')
ylab <- paste( 'PC2 (', pc2, lable,')', sep='')
#par(mar=c(1,1,1,1))
tiff( '"""+output_img+"""',width = 1100, height = 1100,res=150)
plot(d$X1, d$X2, pch=21, bg=sample.bg,  xlab=xlab, ylab=ylab)
abline(h=0,v=0)
text(d$X1, d$X2,d$pc.vector.number,pos=3,cex=0.6,offset=0.5)
dev.off()
		""")
	Rcmd.close()
	subprocess.call( "R CMD BATCH plot_pca.r", stdout=subprocess.PIPE, shell=True )	
	subprocess.call( "rm -rf *.r *.Rout .RData Rplots.pdf", stdout=subprocess.PIPE, shell=True )
	
if __name__ == "__main__": __main__()