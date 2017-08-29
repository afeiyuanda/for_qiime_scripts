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
#sampleID	primerID LibraryID	Time	Locate
ND1	0	Mix1	TimeA	High
ND2	1	Mix1	TimeB	High
ND3	2	Mix1	TimeA	Low
ND4	3	Mix1	TimeB	Low

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
		lib_name = re.split( ',|，', re.split('\s+', line)[2] )[0]
		if lib_name not in lib_dirs:
			lib_dirs.append( lib_name )
			subprocess.call( "mkdir -p " + outdir+'/'+lib_name, stdout=subprocess.PIPE, shell=True )
	return lib_dirs
	
#based on the clo_num column, should be the group name or lib name, starts from 1.
def generate_lines_blocks(lines, clo_num):
	lines_blocks = []
	temp_lines_block = []
	flag_list = []
	for line in lines:
		flag = line.rstrip().split('\t')[int(clo_num)-1]
		if (flag not in flag_list) and len(temp_lines_block)!=0:
			lines_blocks.append( temp_lines_block )
			flag_list.append( flag )
			temp_lines_block = []
			temp_lines_block.append( line )
		else:
			temp_lines_block.append( line )
			if flag not in flag_list:
				flag_list.append( flag )
	lines_blocks.append( temp_lines_block )
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
			flag = re.split('\s+', line)[1]
			sampleid = re.split('\s+', line)[0]
			groups = re.split('\s+', line.strip())[3:]
			for br_pr_line in br_pr_top50:
				br_pr_list = re.split('\s+', br_pr_line.rstrip())
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
			flag = re.split('\s+', line)[1]
			sampleid = re.split('\s+', line)[0]
			groups = re.split('\s+', line.strip())[3:]
			for br_pr_line in br_pr_top50:
				br_pr_list = re.split('\s+', br_pr_line.rstrip())
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
	groups = re.split('\s+', open(sample_primerid_lib).readlines()[0].strip())[3:]
	title_list = ['#SampleID','BarcodeSequence','LinkerPrimerSequence','ReversePrimer']
	title_list.extend(groups)
	title_list.append('Description')
	title_line = '\t'.join(title_list) + '\n'
	big_lib_name_list = []
	for block in lines_blocks:
		lib_name_list = re.split( ',|，', re.split('\s+', block[0])[2] )
		lib_name = lib_name_list[0]
		mapping_file = seqs_dir + '/' + lib_name + '_mapping.txt'
		mapping_file1 = seqs_dir + '/' + lib_name + '_mapping_good.txt'
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
		if lib_name_list not in big_lib_name_list:
			big_lib_name_list.append(lib_name_list)
	return big_lib_name_list

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
	barcode_len = len(re.split('\s+', lines[1])[1])
	#primer_flag==0:Contains Paired Primer
	#primer_flag==2:No Primer
	#primer_flag==1:Single Primer
	if len(re.split('\s+',lines[1].strip())) == 2:
		primer_flag = 2
	elif len(re.split('\s+',lines[1].strip())) == 4:
		primer_flag = 0
	else:
		primer_flag = 1
	
	generate_lib_dirs(outdir, sample_primer_lib)
	seqs_dir = outdir
	big_lib_name_list = generate_mapping_files( seqs_dir, sample_primer_lib, barcode_primer, primer_flag )
	print big_lib_name_list
	
	for lib_name_list in big_lib_name_list:
		f1 = lib_name_list[0]
		f2 = lib_name_list[1]
		os.system('ln -sf '+indir+'/'+f1+' '+seqs_dir+'/'+f1)
		os.system('ln -sf '+indir+'/'+f2+' '+seqs_dir+'/'+f1)#这里f1是没问题的
		print '-----------------------------------Begin to split %s -----------------------------------' % f1
		mapping_file_path = seqs_dir+'/'+f1+'*.txt '
		subprocess.call( "mv " + mapping_file_path + seqs_dir+'/'+f1 , stdout=subprocess.PIPE, shell=True )
			
		print '      --------------Run Flash --------------' 
		#FLASH cmd: flash Mix1_1.fastq.gz  Mix1_2.fastq.gz  -o mix1 -M 310  -t 1 -x 0.5
		fq = seqs_dir+'/'+f1+'/'+f1+'.extendedFrags.fastq'
		if os.path.exists(fq):
			print "%s already exists! Jump to the next step!" % fq
		else:
			cmd = "flash " + ' '.join([seqs_dir+'/'+f1+'/'+i for i in lib_name_list]) + ' -d ' + seqs_dir+'/'+f1 + ' -o ' + f1 + ' -M 310 -t 4'
			print "[[FLASH CMD]]: %s.\nMerged file is %s." % (cmd, fq)
			subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
		
		print '      --------------Edited by My Scrpit grep_from_miseq_fastq_20150115.py--------------' 
		#python ~/scripts/grep_seqs_by_barcode/grep_from_miseq_fastq_20150115.py lib1.extendedFrags.fastq mapping.txt grep_seqs.fastq
		mapping_file = seqs_dir+'/'+f1+'/'+f1 + '_mapping.txt'
		mapping_file1 = seqs_dir+'/'+f1+'/'+f1 + '_mapping_good.txt'
		
		if int(primer_flag) == 0:
			cmd = 'python ' + sys.path[0] + '/grep_seqs_by_barcode/grep_from_miseq_fastq_20150115.py ' + \
			fq+' '+ mapping_file + ' ' + seqs_dir+'/'+f1
		else:
			cmd = 'python ' + sys.path[0] + '/grep_seqs_by_barcode/grep_from_miseq_fastq_20150115_no_primers.py ' + \
			fq+' '+ mapping_file + ' ' + seqs_dir+'/'+f1
		print "[[GREP CMD]]: %s" % cmd
		subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
			
		print '      --------------Fastq Convert To Fasta --------------' 
		#perl ~/software/NGSQCToolkit_v2.3/Format-converter/FastqToFasta.pl lib1_edited.fq
		subprocess.call( 'perl '+sys.path[0]+'/FastqToFasta.pl -i ' + seqs_dir+'/'+f1+'/greped.fq', stdout=subprocess.PIPE, shell=True )
		
		print '      --------------Run split_libraries.py --------------' 
		#split_libraries.py  -f lib1_edited_1line.fna -m mapping1.txt -b 12 -o split_library_output
		input_fasta = seqs_dir+'/'+f1+'/greped.fq_fasta'
		output_dir = seqs_dir+'/'+f1+'/split_library_output'
		if int(primer_flag) == 0:
			cmd = 'split_libraries.py -l 100 -f '+input_fasta+' -b '+str(barcode_len)+' -m ' + mapping_file1 + ' -o ' + output_dir
		else:
			cmd = 'split_libraries.py -l 100 --disable_primers -f '+input_fasta+' -b '+str(barcode_len)+' -m ' + mapping_file1 + ' -o ' + output_dir
		print "[[split_libraries.py CMD]]: %s" % cmd
		subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
		
		print '      --------------Chimera check and Filter --------------' 
		#identify_chimeric_seqs.py -i split_output/seqs.fna -m usearch61 -o usearch_checked_chimeras -r $reference_seqs
		subprocess.call( 'identify_chimeric_seqs.py -i '+output_dir+'/seqs.fna -m usearch61 -o '+seqs_dir+'/'+f1+'/usearch_checked_chimeras -r $reference_seqs', stdout=subprocess.PIPE, shell=True )
		#filter_fasta.py -f split_output/seqs.fna -o chimera_filtered_seqs.fna -s usearch_checked_chimeras/chimeras.txt -n
		subprocess.call( 'filter_fasta.py -f '+output_dir+'/seqs.fna -o '+'/'+seqs_dir+'/'+f1+'/chimera_filtered_seqs.fna -s '+\
		seqs_dir+'/'+f1+'/usearch_checked_chimeras/chimeras.txt -n', stdout=subprocess.PIPE, shell=True )
		print '-----------------------------------End of split %s -----------------------------------\n' % seqs_dir+'/'+f1
		
	subprocess.call( 'cat '+seqs_dir+'/*/usearch_checked_chimeras/chimeras.txt > '+outdir+'/total_chimeras.txt', stdout=subprocess.PIPE, shell=True )
	subprocess.call( 'cat '+seqs_dir+'/*/chimera_filtered_seqs.fna > '+outdir+'/final_seqs.fna', stdout=subprocess.PIPE, shell=True )
	mapping_file_list = glob.glob(seqs_dir+'/*/*_mapping.txt')
	final_mapping_file = outdir+'/mapping.txt'
	final_mapping_file_handle = open(final_mapping_file, 'w')
	final_mapping_file_handle.write(generate_new_map_line(open(mapping_file_list[0]).readlines()[0]))
	for f in mapping_file_list:
		lines = readfile(f)
		for line in lines:
			final_mapping_file_handle.writelines(generate_new_map_line(line))
	final_mapping_file_handle.close()
	
	subprocess.call( 'cat '+seqs_dir+'/*/grep_seqs_num_log.txt > '+outdir+'/grep_seqs_num_log.txt', stdout=subprocess.PIPE, shell=True )
	print '-----You should find your clean seqs in final_seqs.fna----'
	#os.system('rm ' + outdir+'*.fastq')
	
	for lib_name_list in big_lib_name_list:
		f1 = lib_name_list[0]
		os.system('rm -r ' + seqs_dir+'/'+f1)
	
if __name__ == "__main__": __main__()