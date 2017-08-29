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
			Usage:	python average_collated_alpha_index.py collated_alpha_dir legend_pos
			outputs also located in collated_alpha_dir, but with suffix .averaged
			legend pos:¡±bottom¡±,¡±bottomleft¡±,¡±left¡±,¡±topleft¡±,¡±top¡±,¡±topright¡±,¡±right¡±,¡±bottomright¡±,¡±center¡±
		"""

def lines2array( lines_block ):
	#you'd better don't use the same variable name as in the main function
	#change mulitiple lines into array
	arr = []
	for line in lines_block:
		line_list = line.rstrip().split( '\t' )[3:]
		#line_int_list = []
		#for i in line_list:
		#	line_int_list.append( float(i) )
		#arr.append( line_int_list )
		arr.append( line_list )
	return arr
	
def trans_array( arr ):
	rows = len( arr )
	columns = len( arr[0] )
	trans_arr = [[r[col]for r in arr] for col in xrange(columns)] 
	return trans_arr

def float_list( list ):
	return [float(i) for i in list]
		
#in the collated alpha index files, there are many n/a, so we need do some extra operation
def average( list ):
	if 'n/a' in list:
		return 'NA'
	else:
		flist = float_list( list )
		return float(sum(flist))/len(flist)
		
def average_every_column( lines_block ):
	if len( lines_block ) == 1:
		return lines_block[0]
	else:
		arr = lines2array( lines_block )
		trans_arr = trans_array( arr )
		average_list = [ average(i) for i in trans_arr ]
		average_str_list = [ str(i) for i in average_list ]
		average_line = '\t'.join( average_str_list ) + '\n'
	return average_line	

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
		legend_pos = sys.argv[2]
	except:
		usage()
		sys.exit(1)
	
	collated_files = glob.glob( collated_alpha_dir + '/*.txt' )#the files contain the dir path already.
	for file in collated_files:
		f = open( file )
		lines = f.readlines()
		f.close()
		output_ul = file + '.averaged'
		outf = open( output_ul, 'w' )
		title_line = 'numsampled' + '\t' + '\t'.join( lines[0].split('\t')[3:] )
		outf.writelines( title_line )
		flag = [] #record iteration
		for line in lines[1:]:
			if line.split('\t')[2] not in flag:
				flag.append( line.split('\t')[2] )
			else:
				break
		print flag
		new_lines = lines[1:]
		lines_blocks = [new_lines[i:i+len(flag)] for i in range(0, len(new_lines), len(flag))]
		for block in lines_blocks:
			numsampled = block[0].split('\t')[1]
			average_line = numsampled + '\t' + average_every_column( block )
			outf.writelines( average_line )
		outf.close()
	
	#os.chdir(collated_alpha_dir)
	averaged_files = glob.glob( collated_alpha_dir+'/*.averaged' )
	for file in averaged_files:
		Rcmd = open(collated_alpha_dir+"/plot_alpha_index.r", 'w')
		ylab = get_basename(file).split('.')[0]
		output_img= collated_alpha_dir+'/'+ylab + '.png'
		
		Rcmd.writelines(
		"""
library(RColorBrewer)
mycol <- c( brewer.pal(9,"Set1"),brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(8,"Accent"))
mycol <- rep(mycol, 3)
data<-read.table(file='"""+ file +"""', header=T, row.names=1, check.names=F)
new_data <- read.table(file='"""+ file +"""', header=T, row.names=1)
new_data[is.na(new_data)] <- 0
png('"""+ output_img + """',res = 300, width = 3500, height = 3000)
#tiff('"""+ output_img + """',res = 300, width = 3500, height = 3000, compression = "lzw")
#pdf('"""+ output_img + """')
#par(mar=par()$mar+c(0,0,0,6))
plot(x=row.names(data), y=data[,1], ylim=c(min(data[!is.na(data)]), max(new_data)), xlim=c(min(as.numeric(row.names(data))), max(as.numeric(row.names(data)))*(1+0.02)), xlab="Number of Seqs Sampled",ylab='"""+ ylab +"""', type="n", font.lab=1)
lapply( 1:ncol(data), function(i) lines(x=row.names(data), y=data[,i], col=mycol[i], lwd=2.5) )
num_sample = dim(data)[2]
legend('"""+legend_pos+"""',legend=colnames(data), col=mycol[1:num_sample],lty=1,lwd=2.5,cex=1,ncol=2,horiz=F)
#put legend outside the plot region,NOTICE the par parameter!!!!!
#legend(legend=colnames(data), col=mycol[1:num_sample],lty=1,lwd=2.5,cex=1,x=max( as.numeric( row.names(data) ))*(1+0.05),y=min(data, na.rm =T),yjust=0,xpd=TRUE,bty='n')

dev.off()
		""")
		Rcmd.close()
		subprocess.call( "R CMD BATCH "+ collated_alpha_dir+"/plot_alpha_index.r", stdout=subprocess.PIPE, shell=True )
		#subprocess.call( "rm -rf "+collated_alpha_dir+"/*.r "+collated_alpha_dir+"/*.Rout "+collated_alpha_dir+"/*.RData "+collated_alpha_dir+"/Rplots.pdf ", stdout=subprocess.PIPE, shell=True )
		subprocess.call( "convert " + output_img + ' ' + collated_alpha_dir+'/'+ylab + '.png' , stdout=subprocess.PIPE, shell=True )
	
if __name__ == "__main__": __main__()