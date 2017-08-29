#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.02
Usage:	python plot_venn_diagram.py   sorted_otu_table_by_Treatment.txt
'''

import sys,subprocess,os

def usage():
    print """
			Usage:	python plot_venn_diagram.py   sorted_otu_table_by_Treatment.txt
			OTU table should be Grouped firstly!
			Input file like this:
				# Constructed from biom file
				#OTU ID	M	D	F
				3951715	0.0	0.0	2.0
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
	new_lines.insert(0,lines[1][1:]) #[1:] means ignore the # at the beginning of second line
	return new_lines

#line:3951715	0.0	0.0	2.0,the first col is OTU ID
#Then return a new line
def map_otuid_to_samples( line ):
	line_list = line.rstrip().split('\t')
	otuid = line_list[0]
	reads_num_list = line_list[1:]
	new_num_list = []
	for i in reads_num_list:
		if i == '0.0':
			new_num_list.append('NA')
		else:
			new_num_list.append(otuid)
	return '\t'.join( new_num_list ) + '\n'

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
		grouped_otu_table = sys.argv[1]
	except:
		usage()
		sys.exit(1)
	
	otu_table_dir = get_dirname(grouped_otu_table)
	otu_table_name = get_basename(grouped_otu_table)
	
	output = otu_table_dir + '/' + otu_table_name + '.formated'
	lines = readfile( grouped_otu_table )
	sample_list = lines[0].rstrip().split('\t')[1:]
	sample_num = len( sample_list )
	sample_name_line = '\t'.join( sample_list ) + '\n'
	otu_lines = lines[1:]
	
	outf = open( output, 'w' )
	outf.writelines( sample_name_line )
	for line in otu_lines:
		new_line = map_otuid_to_samples( line )
		outf.writelines( new_line )
	outf.close()
	
	rscript_name = otu_table_dir+"/plot_venn_diagram_" + str(sample_num) + '.r'
	venn_name = otu_table_dir+'/'+otu_table_name.split('.')[0] + '_' + str(sample_num) + '.png'
	Rcmd = open(rscript_name, 'w')
	if sample_num == 5:
		Rcmd.writelines(
		"""
library(VennDiagram)
file <- '"""+ output +"""'
read.table(file,header=TRUE,sep="\\t",check.names=F) -> venndata
dataf<-as.list(venndata)
new_dataf <- lapply(dataf, function(x) x[!is.na(x)])
mycol = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
venn.plot <- venn.diagram(
new_dataf,
imagetype = "png",
filename = '"""+ venn_name +"""',
col = "black",
fill = mycol,
alpha = 0.50,
#cex = c(2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 0.8, 1.5, 0.8, 1.5, 0.8, 1.5, 0.8,1.5, 0.8, 1.5, 0.55, 1.5, 0.55, 1.5, 0.55, 1.5, 0.55, 1.5, 0.55, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5),
cat.col = mycol,
cat.cex = 1,
cat.fontface = "bold",
margin = 0.05)
		""")
	else:
		Rcmd.writelines(
		"""
library(VennDiagram)
file <- '"""+ output +"""'
read.table(file,header=TRUE,sep="\\t",check.names=F) -> venndata
dataf<-as.list(venndata)
new_dataf <- lapply(dataf, function(x) x[!is.na(x)])
mycol = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
venn.plot <- venn.diagram(
new_dataf,
imagetype = "png",
filename = '"""+ venn_name +"""',
col = "black",
fill = mycol[1:"""+str(sample_num)+"""],
alpha = 0.50,
#cex = c(2.5, 2.5, 2.5, 2.5, 2.5, 2, 0.8, 2, 0.8, 2, 0.8, 2, 0.8,2, 0.8, 2, 0.55, 2, 0.55, 2, 0.55, 2, 0.55, 2, 0.55, 2, 2, 2, 2, 2, 2.5),
cat.col = mycol[1:"""+str(sample_num)+"""],
cat.cex = 1,
cat.fontface = "bold",
margin = 0.05)
		""")
		
	Rcmd.close()
	subprocess.call( "R CMD BATCH " + rscript_name, stdout=subprocess.PIPE, shell=True )
	subprocess.call( "convert " + venn_name + ' '+otu_table_dir+ '/venn_diagram_' + otu_table_name + '_' + str(sample_num) + '.png', stdout=subprocess.PIPE, shell=True )
	#subprocess.call( "rm -rf *.r *.Rout .RData Rplots.pdf", stdout=subprocess.PIPE, shell=True )
	
if __name__ == "__main__": __main__()