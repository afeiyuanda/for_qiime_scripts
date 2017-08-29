#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.11.19, 2014.12.04,
		2014.12.10:can batch process
		2015.02.02:if only two groups,we shoud define two color directly!
Usage:	python barplot_average_taxo_abundance.py mapping_file  original_taxo_file  output_prefix
'''

import sys,copy,subprocess,glob

def usage():
    print """
			Usage:	python barplot_average_taxo_abundance.py wf_taxa_summary/ mapping_file.txt group_flag
			Your should notice that different group have different NO. of samples!
		"""
		
#This fun return a list,list[0] is groups and its samples formed dic, list[1] is groups empty dic, list[2] is the min sample NO. in all groups
def get_groups_info(mapping_file_path, group_flag):
	f = open( mapping_file_path )
	lines = f.readlines()
	groups_dic = {}
	groups = []
	group_flag_id = lines[0].split('\t').index(group_flag)
	for line in lines:
		if line.startswith('#'):
			continue
		else:
			group = line.rstrip().split('\t')[group_flag_id]
			if group not in groups:
				groups.append( group )
	#print groups	
	for group in groups:
		groups_dic[group] = []
	empty_groups_dic = copy.deepcopy(groups_dic) #we use copy.deepcopy, intend to make a new dic that will not change when groups_dic change
	
	for group in groups:
		for line in lines:
			if line.rstrip().split('\t')[group_flag_id]==group:
				groups_dic[group].append( line.rstrip().split('\t')[0] )
	for key in groups_dic.keys():
		if len(groups_dic[key]) == 1:
			del groups_dic[key]
	print groups_dic		
	min_samples = 10000 #keep the min sample NO. in all groups
	max_samples = 0
	for key in groups_dic.keys():
		if len( groups_dic[key] ) < min_samples:
			min_samples = len( groups_dic[key] )
		if len( groups_dic[key] ) > max_samples:
			max_samples = len( groups_dic[key] )
			
	result_list = []
	result_list.append( groups_dic )
	result_list.append( empty_groups_dic )
	result_list.append( min_samples )
	result_list.append( max_samples )
	return result_list
		
def __main__():

	try:
		wf_taxa_summary = sys.argv[1]
		mapping_file = sys.argv[2]
		group_flag = sys.argv[3]
	except:
		usage()
		sys.exit(1)
	groups_dic = get_groups_info(mapping_file,group_flag)[0]
	print groups_dic
	groups_no = len(groups_dic.keys())
	#{'M': ['M1', 'M2', 'M3'], 'D': ['D1', 'D2', 'D3'], 'F': ['F1', 'F2', 'F3']}
	empty_groups_dic = get_groups_info(mapping_file,group_flag)[1]
	#{'M': [], 'D': [], 'F': []}
	min_sample = get_groups_info(mapping_file,group_flag)[2]
	
	level_files = glob.glob( wf_taxa_summary + '/L*')
	print level_files
	for level_file in level_files:
		original_taxo_file = level_file + '/average.xls'
		f = open( original_taxo_file )
		lines = f.readlines()
		samples = lines[0].rstrip().split('\t')
		for group in groups_dic.keys():
			for sample in groups_dic[group]:
				empty_groups_dic[group].append( samples.index( sample ) )
		#empty_groups_dic now changed as following:
		#{'M': [7, 8, 9], 'D': [1, 2, 3], 'F': [4, 5, 6]}
		output_prefix = 'taxo'
		reshap_taxo_file = open( output_prefix + '.txt', 'w' )
		title_line = "\t".join(["mygroups", "taxo"] + ["sample"+str(i) for i in range(min_sample)])
		reshap_taxo_file.writelines( title_line + '\n')
		for group in groups_dic.keys():
			for line in lines[1:]:
					new_line = group + '\t' + line.split('\t')[0] + '\t' + '\t'.join( [line.rstrip().split('\t')[i] for i in empty_groups_dic[group][0:min_sample] ] ) + '\n'
					reshap_taxo_file.writelines( new_line )
		reshap_taxo_file.close()
		
		barplot_img = level_file + '/taxo_relative_abundance.png'
		if min_sample == 1:
			print "In every group only have one sample!"
			Rcmd = open(level_file+"/plot_bargraphCI.r", 'w')
			Rcmd.writelines(
			"""
library(RColorBrewer)
png(file = '"""+ barplot_img +"""',width=1000, height=0,res=150)
#par(xpd=T)
par(xpd=T,mar=c(6,3,0,0), oma=c(1,1,0,0))
taxo_data <- read.table('"""+ original_taxo_file +"""', header=T, sep="\t",row.names=1) 
mp <- barplot(t(taxo_data), beside=TRUE, col=brewer.pal(""" + str(groups_no) + ""","Set1"), ylim = c(0, max( taxo_data )+0.1),lwd = 1, names.arg=NULL,axisnames=FALSE) 
x <- apply(mp,2,max)
text(cex=0.8,x, y=-0.01, lab=row.names(taxo_data), srt=45, adj=1)
mtext(1, text ="Taxonomy", line = 6, cex=1) #cex control the font size, line control the distance between text and the x axis,negtive=up, positive=down
mtext("Relative abundance +/- SD(%)",side=2,line=2,font=0.3,cex=1)  #add y lable
legend("topright", names(taxo_data), bty="n",horiz = T, fill=brewer.pal(""" + str(groups_no) + ""","Set1"))
dev.off()
			"""
							)
			Rcmd.close()
		
		elif groups_no == 2:
			print "There are 2 groups!"
			Rcmd = open(level_file+"/plot_bargraphCI.r", 'w')
			Rcmd.writelines(
				"""
library(reshape)
library(sciplot)
library(RColorBrewer)
png(file = '"""+ barplot_img +"""',width=1000, height=1000,res=150)
#tiff(file = '"""+ barplot_img +"""',width=1000, height=1000,res=150)
#pdf(file = '"""+ barplot_img +"""')
oldmar <- par()$mar
par(xpd=T,mar=c(6,3,0,0), oma=c(1,1,0,0))
#par(xpd=T,mar=c(6,6,0,0), oma=c(5,5,3,3))
#par(xpd=T)
d <- read.table("taxo.txt",header=T, sep='\\t')
d$taxo <- as.factor(d$taxo)
newtest <- melt(d)
two_col <- c("#E41A1C", "#377EB8")
mp <- bargraph.CI(taxo, value, group = mygroups,err.width=if(length(levels(newtest$taxo))>10) 0.02 else .05, data = newtest, col=two_col, err.col = "black", ci.fun = function(x) c(mean(x)-sd(x), mean(x)+sd(x)), ylim = c(0, max( newtest$value )+0.1),lwd = 1, names.arg=NULL,axisnames=FALSE,lc=F)
x <- apply(mp$xvals,2,max)
text(cex=0.8,x, y=-0.01, lab=levels(newtest$taxo), srt=45, adj=1)
mtext(1, text ="Taxonomy", line = 6, cex=1) #cex control the font size, line control the distance between text and the x axis,negtive=up, positive=down
mtext("Relative abundance +/- SD(%)",side=2,line=2,font=0.3,cex=1)  #add y lable
legend("topright", legend = levels(newtest$mygroups), bty = "n", horiz = T, fill = two_col)
dev.off()
par(oldmar)
			"""
						)
			Rcmd.close()
		
		else:
			print "There are more than 2 groups!" 
			Rcmd = open(level_file+"/plot_bargraphCI.r", 'w')
			Rcmd.writelines(
				"""
library(reshape)
library(sciplot)
library(RColorBrewer)
png(file = '"""+ barplot_img +"""',width=1000, height=1000,res=150)
#tiff(file = '"""+ barplot_img +"""',width=1500, height=1500,res=150)
#pdf(file = '"""+ barplot_img +"""')
par(xpd=T,mar=c(6,3,0,0), oma=c(1,1,0,0))
#par(xpd=T,mar=c(6,6,0,0), oma=c(5,5,3,3))
#par(xpd=T)
d <- read.table("taxo.txt",header=T, sep='\\t')
d$taxo <- as.factor(d$taxo)
newtest <- melt(d)
mp <- bargraph.CI(taxo, value, group = mygroups,err.width=if(length(levels(newtest$taxo))>10) 0.02 else .05, data = newtest, col=brewer.pal(""" + str(groups_no) + ""","Set1"), err.col = "black", ci.fun = function(x) c(mean(x)-sd(x), mean(x)+sd(x)), ylim = c(0, max( newtest$value )+0.1),lwd = 1, names.arg=NULL,axisnames=FALSE,lc=F)
x <- apply(mp$xvals,2,max)
text(cex=0.8,x, y=-0.01, lab=levels(newtest$taxo), srt=45, adj=1)
mtext(1, text ="Taxonomy", line = 6, cex=1) #cex control the font size, line control the distance between text and the x axis,negtive=up, positive=down
mtext("Relative abundance +/- SD(%)",side=2,line=2,font=0.3,cex=1)  #add y lable
legend("topright", legend = levels(newtest$mygroups), bty = "n", horiz = T, fill = brewer.pal(""" + str(groups_no) + ""","Set1"))
dev.off()
			"""
							)
	
			Rcmd.close()
		
		subprocess.call( "R CMD BATCH "+level_file+"/plot_bargraphCI.r", stdout=subprocess.PIPE, shell=True )
		subprocess.call( "rm -rf *.r *.Rout .RData Rplots.pdf", stdout=subprocess.PIPE, shell=True )
	
	
if __name__ == "__main__": __main__()