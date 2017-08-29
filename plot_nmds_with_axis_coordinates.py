#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.23
Usage:	python plot_nmds_with_axis_coordinates.py nmds_weighted_unifrac_otu_table_even1000.txt mapping.txt Treatment legend_pos
'''

import sys,glob,subprocess,os

def usage():
    print """
			Usage:	python plot_nmds_with_axis_coordinates.py nmds_weighted_unifrac_otu_table_even1000.txt mapping.txt Treatment legend_pos
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
	Rcmd = open("plot_nmds.r", 'w')	
	Rcmd.writelines(
		"""
library(RColorBrewer)
col_list <- c( brewer.pal(9,"Set1"),brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(8,"Accent"))
new_pch <- c(c(21:22),c(24:25))
pch_list <- rep( new_pch,5 )

d <- read.table('"""+ output_name +"""',header=T,sep='\\t')
samplelist = d$group # "F" "F" "F" "M" "M" "M" "D" "D" "D"
sample2factor_num = as.numeric( samplelist ) #2 2 2 3 3 3 1 1 1

sample.bg=col_list[sample2factor_num]
legend_annot = levels(d$group)

legend_factor = as.numeric(levels( as.factor(sample2factor_num )))
lengend.bg=col_list[legend_factor]
xlab <- 'NMDS1'
ylab <- 'NMDS2'
#par(mar=c(1,1,1,1))
tiff( '"""+output_img+"""',width = 1100, height = 1100,res=150)
#plot(d$NMDS1, d$NMDS2, pch=pch_list[sample2factor_num], bg=sample.bg,  xlab=xlab, ylab=ylab,xlim=c(min(d$X1)0.05, max(d$X2)+0.05))
plot(d$NMDS1, d$NMDS2, pch=pch_list[sample2factor_num], bg=sample.bg,  xlab=xlab, ylab=ylab)
abline(h=0,v=0)
legend('"""+ lengend_pos +"""',as.graphicsAnnot( legend_annot ), pch=pch_list[legend_factor],pt.bg=lengend.bg,cex=0.8,pt.cex=1.1)
text(d$NMDS1, d$NMDS2,d$samples,pos=3,cex=0.6,offset=0.5)
dev.off()
		""")
	Rcmd.close()
	subprocess.call( "R CMD BATCH plot_nmds.r", stdout=subprocess.PIPE, shell=True )
	subprocess.call( "convert " + output_img + ' ' + pcoa_file.split('.')[0] + '.png', stdout=subprocess.PIPE, shell=True )
	subprocess.call( "rm -rf *.r *.Rout .RData Rplots.pdf", stdout=subprocess.PIPE, shell=True )	
	
if __name__ == "__main__": __main__()