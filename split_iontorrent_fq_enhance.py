#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.15
Usage:	python split_librarys_batch.py top50_barcodes_primers.txt sample_primerid_lib.txt
'''

import sys,os,glob,subprocess,re,argparse

def usage():
    print """
Usage:	python """+sys.argv[0]+ """ 
-indir raw_fq_dir 
-barcode_primer /path/to/barcode_primer_file.txt 
-sample_primer_lib 	/path/to/sample_primer_lib.txt 
-outdir /path/to/outdir

datas name like:Mix10_1.fastq.gz  Mix12_1.fastq.gz

---------------top50_barcodes_primers.txt----------------
#PrimerID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer
0	TCCCTTGTCTCC	GTGCCAGCMGCCGCGGTAA	GGACTACHVGGGTWTCTAAT
1	ACGAGACTGATT	GTGCCAGCMGCCGCGGTAA	GGACTACHVGGGTWTCTAAT

---------------sample_primerid_lib.txt----------------
#sampleID	primerID fastq_name	Time	Locate
ND1	0	1106yx.fastq	TimeA	High
ND2	1	1106yx.fastq	TimeB	High
ND3	2	1106yx.fastq	TimeA	Low
ND4	3	1106yx.fastq	TimeB	Low

#primer_flag==0:Contains Paired Primer
#primer_flag==2:No Primer
#primer_flag==1:Single Primer
		"""

def readfile( file_ul ):
	f = open( file_ul )
	lines = f.readlines()
	new_lines = []
	for line in lines:
		if line.startswith('#'):
			continue
		else:
			new_lines.append( line )
	return new_lines
	
def generate_lib_dirs(outdir, sample_primerid_lib):
	f = open(sample_primerid_lib, 'r')
	lines = f.readlines()
	f.close()
	lib_dirs = []
	for line in lines:
		if line.startswith('#'):
			continue
		lib_name = re.split('\s+', line)[2]
		if lib_name not in lib_dirs:
			lib_dirs.append( lib_name )
			subprocess.call( "mkdir -p " + outdir+'/'+lib_name, stdout=subprocess.PIPE, shell=True )
	return lib_dirs

#based on the clo_num column, should be the group name or lib name, starts from 1.
#根据Group name或者lib name对样品编号文件进行分块，注意：col num是从1开始的
def generate_lines_blocks(lines, clo_num):
	lines_blocks = []
	temp_lines_block = []
	flag_list = []
	for line in lines:
		flag = line.rstrip().split('\t')[int(clo_num)-1]
		if flag not in flag_list:
			flag_list.append(flag)
	for flag in flag_list:
		for line in lines:
			line_flag = line.rstrip().split('\t')[int(clo_num)-1]
			if flag == line_flag:
				temp_lines_block.append(line)
		lines_blocks.append(temp_lines_block)
		temp_lines_block = []
	#print flag_list
	return lines_blocks

#primer_flag==0:Contains Paired Primer
#primer_flag==2:No Primer
#primer_flag==1:Single Primer
def generate_mapping_file_lines( sa_pr_lib_lines, br_pr_top50, primer_flag ):
	mapping_file_lines = []
	for line in sa_pr_lib_lines:
		if line.startswith('#'):
			pass
		else:
			flag = line.split('\t')[1]
			sampleid = line.split('\t')[0]
			groups = line.strip().split('\t')[3:]
			for br_pr_line in br_pr_top50:
				br_pr_list = br_pr_line.rstrip().split('\t')
				primer_id = br_pr_list[0]
				if flag == primer_id and int(primer_flag) == 0:
					print 'primer_flag == 0'
					mapping_file_line = sampleid + '\t' + '\t'.join( br_pr_list[1:] ) + '\t' + '\t'.join(groups)+'\t' + primer_id + '\n'
					break
				elif flag == primer_id and int(primer_flag) == 1:
					print 'primer_flag == 1'
					mapping_file_line = sampleid + '\t' + '\t'.join( br_pr_list[1:] ) + '\t\t' + '\t'.join(groups)+'\t' + primer_id + '\n'
					break
				elif flag == primer_id and int(primer_flag) == 2:
					print 'primer_flag == 2'
					mapping_file_line = sampleid + '\t' + '\t'.join( br_pr_list[1:] ) + '\t\t\t' + '\t'.join(groups)+'\t' + primer_id + '\n'
					break
			mapping_file_lines.append( mapping_file_line )
	return mapping_file_lines

def generate_mapping_file_lines1( sa_pr_lib_lines, br_pr_top50, primer_flag ):
	mapping_file_lines = []
	for line in sa_pr_lib_lines:
		if line.startswith('#'):
			pass
		else:
			flag = line.split('\t')[1]
			sampleid = line.split('\t')[0]
			groups = line.strip().split('\t')[3:]
			for br_pr_line in br_pr_top50:
				br_pr_list = br_pr_line.rstrip().split('\t')
				primer_id = br_pr_list[0]
				if flag == primer_id and int(primer_flag) == 0:
					print 'primer_flag == 0'
					mapping_file_line = sampleid + '\t' + '\t'.join( br_pr_list[1:] ) + '\t' + '\t'.join(groups)+'\t' + primer_id + '\n'
					break
				elif flag == primer_id and int(primer_flag) == 1:
					print 'primer_flag == 1'
					mapping_file_line = sampleid + '\t' + '\t'.join( br_pr_list[1:] ) + '\tGGACTACHVGGGTWTCTAAT\t' + '\t'.join(groups)+'\t' + primer_id + '\n'
					break
				elif flag == primer_id and int(primer_flag) == 2:
					print 'primer_flag == 2'
					mapping_file_line = sampleid + '\t' + '\t'.join( br_pr_list[1:] ) + '\tGTGCCAGCMGCCGCGGTAA\tGGACTACHVGGGTWTCTAAT\t' + '\t'.join(groups)+'\t' + primer_id + '\n'
					break
			mapping_file_lines.append( mapping_file_line )
	return mapping_file_lines
	
def generate_mapping_files( seqs_dir, sample_primerid_lib, top50_barcodes_primers, primer_flag ):
	sa_pr_lib_lines = readfile( sample_primerid_lib )
	br_pr_top50 = readfile( top50_barcodes_primers )
	lines_blocks = generate_lines_blocks( sa_pr_lib_lines, 3 )
	groups = open(sample_primerid_lib).readlines()[0].strip().split('\t')[3:]
	title_list = ['#SampleID','BarcodeSequence','LinkerPrimerSequence','ReversePrimer']
	title_list.extend(groups)
	title_list.append('Description')
	title_line = '\t'.join(title_list) + '\n'
	libid_list = []
	for block in lines_blocks:
		libid = block[0].rstrip().split('\t')[2]
		mapping_file = seqs_dir + '/' + libid + '_mapping.txt'
		mapping_file1 = seqs_dir + '/' + libid + '_mapping_good.txt'
		moutput = open( mapping_file, 'w' )
		moutput1 = open( mapping_file1, 'w' )
		moutput.writelines( title_line )
		moutput1.writelines( title_line )
		mapping_file_lines = generate_mapping_file_lines( block, br_pr_top50, primer_flag )
		mapping_file_lines1 = generate_mapping_file_lines1( block, br_pr_top50, primer_flag )
		for mapping_file_line in mapping_file_lines:
			moutput.writelines( mapping_file_line )
		moutput.close()
		for mapping_file_line in mapping_file_lines1:
			moutput1.writelines( mapping_file_line )
		moutput1.close()
		if libid not in libid_list:
			libid_list.append(libid)
	return libid_list

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

#只取mapping文件中第一列，分组列，及最后一列
def generate_new_map_line(line):
	line_list = line.strip().split('\t')
	new_line = line_list[0]+'\t'+'\t'.join(line_list[4:])+'\n'
	return new_line
		
def __main__():
	
	#quick_usage= 'python '+sys.argv[0]+ ' -indir raw_fq_dir -barcode_primer /path/to/barcode_primer_file.txt -sample_primer_lib \
	#/path/to/sample_primer_lib.txt -outdir /path/to/outdir'

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-indir", dest="indir", help="Directory contains all rawdata", required=True);
	newParser.add_argument( "-barcode_primer", dest="barcode_primer", help="File contains all the barcodes and primers seqs", required=True);
	newParser.add_argument( "-sample_primer_lib", dest="sample_primer_lib", help="File contains SampleID and corresponding PrimerID and LibID" );
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	indir = argsDict['indir']
	barcode_primer = argsDict['barcode_primer']
	sample_primer_lib = argsDict['sample_primer_lib']
	outdir = argsDict['outdir']
	subprocess.call( 'dos2unix '+barcode_primer, stdout=subprocess.PIPE, shell=True )
	subprocess.call( 'dos2unix '+sample_primer_lib, stdout=subprocess.PIPE, shell=True )
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	f = open(barcode_primer)
	lines = f.readlines()
	f.close()
	barcode_len = len(lines[1].split('\t')[1])
	#primer_flag==0:Contains Paired Primer
	#primer_flag==2:No Primer
	#primer_flag==1:Single Primer
	if len(re.split('\s+',lines[1].strip())) == 2:
		primer_flag = 2
	elif len(re.split('\s+',lines[1].strip())) == 4:
		primer_flag = 0
	else:
		primer_flag = 1
	
	seqs_dir = outdir
	generate_mapping_files( seqs_dir, sample_primer_lib, barcode_primer, primer_flag )	
	lib_dirs = generate_lib_dirs(outdir, sample_primer_lib)
	print lib_dirs
	
	for lib_dir in lib_dirs:
		print '-----------------------------------Begin to split %s -----------------------------------' % lib_dir 
		os.system('ln -sf '+indir+'/'+lib_dir+' '+seqs_dir+'/'+lib_dir)
		mapping_file_path = seqs_dir+'/'+lib_dir + '*.txt '
		subprocess.call( "mv " + mapping_file_path + ' ' + seqs_dir+'/'+lib_dir , stdout=subprocess.PIPE, shell=True )
	
		if lib_dir.endswith('.gz'):
			os.system('tar -zxvf ' + seqs_dir+'/'+lib_dir+'/'+lib_dir)
			fastq_file = seqs_dir+'/'+lib_dir+'/'+lib_dir.split('.gz')[0]
		else:
			fastq_file = seqs_dir+'/'+lib_dir+'/'+lib_dir
		subprocess.call( "python " + sys.path[0] + '/split_iontorrent_fq.py -fastq '+fastq_file+' -mapping '+seqs_dir+'/'+lib_dir+'/'+lib_dir+'_mapping.txt -outdir '+seqs_dir+'/'+lib_dir, stdout=subprocess.PIPE, shell=True )
			
		print '-----------------------------------End of split %s --------------------------------------\n' % lib_dir
		
	#os.system('mkdir -p '+outdir+'/Combined_Seqs')
	subprocess.call( 'cat '+outdir+'/*/final_seqs.fna > '+outdir+'/final_seqs.fna', stdout=subprocess.PIPE, shell=True )
	#将所有的mapping文件整合起来
	mapping_file_list = glob.glob(seqs_dir+'/*/*_mapping.txt')
	final_mapping_file = outdir+'/mapping.txt'
	final_mapping_file_handle = open(final_mapping_file, 'w')
	final_mapping_file_handle.write(generate_new_map_line(open(mapping_file_list[0]).readlines()[0]))
	for f in mapping_file_list:
		lines = readfile(f)
		for line in lines:
			final_mapping_file_handle.writelines(generate_new_map_line(line))
	final_mapping_file_handle.close()
	
	subprocess.call( 'cat '+seqs_dir+'/*/grep_log.txt > '+outdir+'/grep_log.txt', stdout=subprocess.PIPE, shell=True )
	subprocess.call( 'cat '+seqs_dir+'/*/split_library_log.txt > '+outdir+'/split_library_log.txt', stdout=subprocess.PIPE, shell=True )
	for lib_dir in lib_dirs:
		os.system('rm -rf '+seqs_dir+'/'+lib_dir)
if __name__ == "__main__": __main__()