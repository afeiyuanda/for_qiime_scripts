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
-sample_lib 	/path/to/sample_lib.txt 
-outdir /path/to/outdir

datas name like:1_1.fq.gz  1_2.fq.gz

---------------sample_lib.txt----------------
#SampleID SeqsNames	Time	Locate
ND1	1_1.fq.gz,1_2.fq.gz	TimeA	High
ND2	2_1.fq.gz,2_2.fq.gz	TimeB	High
ND3	3_1.fq.gz,3_2.fq.gz	TimeA	Low
ND4	4_1.fq.gz,4_2.fq.gz	TimeB	Low
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

#只取mapping文件中第一列，分组列，及最后一列
def generate_new_map_line(line):
	line_list = line.strip().split('\t')
	new_line = line_list[0]+'\t'+'\t'.join(line_list[3:])+'\n'
	return new_line
	
def __main__():
	
	#quick_usage= 'python '+sys.argv[0]+ ' -indir raw_fq_dir -sample_lib /path/to/sample_lib.txt -outdir /path/to/outdir'

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-indir", dest="indir", help="Directory contains all rawdata.", required=True);
	newParser.add_argument( "-sample_lib", dest="sample_lib", help="File contains SampleID and corresponding LibID and Group infor." );
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir.", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	indir = argsDict['indir']
	sample_lib = argsDict['sample_lib']
	outdir = argsDict['outdir']
	subprocess.call( 'dos2unix '+sample_lib, stdout=subprocess.PIPE, shell=True )
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	os.system('mkdir '+outdir+'/clean_data')
	sample_lib_lines = open( sample_lib ).readlines()[1:]
	title_line = open( sample_lib ).readlines()[0]
	group_names = re.split('\s+', title_line.strip())[2:]
	for sample_lib in sample_lib_lines:
		sample_lib_list = re.split('\s+', sample_lib.strip())
		sample_id = sample_lib_list[0]
		fq1 = re.split(',|，', sample_lib_list[1])[0]
		fq2 = re.split(',|，', sample_lib_list[1])[1]
		groups = sample_lib_list[2:]
		fq1_outdir = outdir+'/'+fq1+'.dir'
		os.system('mkdir '+fq1_outdir)
		os.system('ln -sf '+indir+'/'+fq1+' '+fq1_outdir)
		os.system('ln -sf '+indir+'/'+fq2+' '+fq1_outdir)
		#mapping.txt的格式如下
		#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	Parallel	Description
		#ML1			ML	L	10714L
		#MS1			MS	S	10714S
		f = open(fq1_outdir+'/mapping.txt', 'w')
		f.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\t'+'\t'.join(group_names)+'\tDescription\n')
		f.write(sample_id+'\t\t\t'+'\t'.join(groups)+'\t'+sample_lib_list[1]+'\n')
		f.close()
		print '-----------------------------------Begin to split %s -----------------------------------' % fq1
		print '      --------------Run Flash --------------' 
		#FLASH cmd: flash Mix1_1.fastq.gz  Mix1_2.fastq.gz  -o mix1 -M 310  -t 1 -x 0.5
		cmd = 'flash ' + fq1_outdir+'/'+fq1+' '+fq1_outdir+'/'+fq2+' -d '+fq1_outdir+' -o ' + fq1 + ' -M 310 -t 4'
		print "[[FLASH CMD]]: %s.\nMerged file prefix is %s." % (cmd, fq1)
		subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
		
		print '      --------------Run split_libraries_fastq.py --------------' 
		#split_libraries_fastq.py -i MIX1.extendedFrags.fastq -m mapping1.txt --sample_ids D1.C1 --barcode_type 'not-barcoded' -o split_library_output -q 30  -p 0.80
		cmd = 'split_libraries_fastq.py -i '+ fq1_outdir+'/'+fq1+'.extendedFrags.fastq --barcode_type not-barcoded -q 20  -p 0.80 --store_demultiplexed_fastq -m ' + fq1_outdir+'/mapping.txt -o ' + fq1_outdir+'/split_library_output --sample_ids ' + sample_id
		print '[[split_libraries_fastq CMD]]:%s' % cmd
		subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
		subprocess.call( 'tail -5 ' + fq1_outdir+'/split_library_output/split_library_log.txt | head -1 >> '+outdir+'/split_seqs_log.txt', stdout=subprocess.PIPE, shell=True )
		os.system('cp '+fq1_outdir+'/split_library_output/seqs.fastq '+outdir+'/clean_data/'+fq1.split('.')[0]+'.fastq')
		
		print '      --------------Chimera check and Filter --------------' 
		#identify_chimeric_seqs.py -i split_library_output/seqs.fna -m usearch61 -o usearch_checked_chimeras -r $reference_seqs
		cmd = 'identify_chimeric_seqs.py -i '+fq1_outdir+'/split_library_output/seqs.fna -m usearch61 -o '+fq1_outdir+'/usearch_checked_chimeras -r /share/nas2/genome/biosoft/QIIME/gg_otus-13_8-release/rep_set/97_otus.fasta'
		print '[[identify_chimeric_seqs CMD]]:%s' % cmd
		subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
		
		#filter_fasta.py -f split_library_output/seqs.fna -o chimera_filtered_seqs.fna -s usearch_checked_chimeras/chimeras.txt -n
		cmd = 'filter_fasta.py -f '+fq1_outdir+'/split_library_output/seqs.fna -o '+fq1_outdir+'/chimera_filtered_seqs.fna -s '+fq1_outdir+'/usearch_checked_chimeras/chimeras.txt -n'
		print '[[filter_fasta CMD]]:%s' % cmd
		subprocess.call( cmd, stdout=subprocess.PIPE, shell=True )
		print '-----------------------------------End of split %s -----------------------------------\n' % fq1_outdir
		
	subprocess.call( 'cat '+outdir+'/*/chimera_filtered_seqs.fna > '+outdir+'/final_seqs.fna', stdout=subprocess.PIPE, shell=True )
	mapping_file_list = glob.glob(outdir+'/*/mapping.txt')
	final_mapping_file = outdir+'/mapping.txt'
	final_mapping_file_handle = open(final_mapping_file, 'w')
	final_mapping_file_handle.write(generate_new_map_line(open(mapping_file_list[0]).readlines()[0]))
	for f in mapping_file_list:
		lines = readfile(f)
		for line in lines:
			final_mapping_file_handle.writelines(generate_new_map_line(line))
	final_mapping_file_handle.close()
	print '-----You should find your clean seqs in final_seqs.fna----'
	#os.system('rm -r ' + outdir+'/*/')
	
if __name__ == "__main__": __main__()
