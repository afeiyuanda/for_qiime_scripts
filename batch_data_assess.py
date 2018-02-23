#!/usr/bin/python  
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.09.30
Usage:	
'''

import os,sys,shutil,glob,re,datetime
import argparse,subprocess

def read_mapping_file(mapping_file_url):
	f = open(mapping_file_url)
	lines = f.readlines()
	f.close()
	d = {}
	for line in lines[1:]:
		sample_id = line.strip().split('\t')[0]
		sample_data = line.strip().split('\t')[-1]
		d[sample_id] = sample_data
	return d

def generate_data_assess_cfg(rawdata_dir, mapping_file_url):
	data_dic = read_mapping_file(mapping_file_url)
	cfg = rawdata_dir + '/data_assess.cfg'
	f = open(cfg, 'w')
	for sample in sorted(data_dic.keys()):
		f.write('Sample\t' + sample + '\n')
		seqs = data_dic[sample]
		fq1 = re.split(',|，',data_dic[sample])[0].split('.gz')[0]
		fq2 = re.split(',|，',data_dic[sample])[1].split('.gz')[0]
		f.write('fq1\t' + rawdata_dir + '/' +fq1+'\n')
		f.write('raw_fq1\t' + rawdata_dir + '/' + fq1 + '\n')
		f.write('fq2\t' + rawdata_dir + '/' + fq2 + '\n')
		f.write('raw_fq2\t' + rawdata_dir + '/' + fq2 + '\n')
	f.close()
	print "Data Assess config file generate successfully at %s !" % (rawdata_dir + '/data_assess.cfg')

def generate_data_assess_cfg1(rawdata_dir, mapping_file_url):
	data_dic = read_mapping_file(mapping_file_url)
	cfg = rawdata_dir + '/data_assess.cfg'
	f = open(cfg, 'w')
	for sample in sorted(data_dic.keys()):
		seqs = data_dic[sample]
		fq1 = re.split(',|，',data_dic[sample])[0].split('.fastq')[0]+'.fastq'
		#fq1 = re.split(',|，',data_dic[sample])[0].split('.gz')[0]
		f.write('Sample\t' + sample + '\n')
		f.write('fq1\t' + rawdata_dir + '/' +fq1+'\n')
		f.write('raw_fq1\t' + rawdata_dir + '/' + fq1 + '\n')
	f.close()
	print "Data Assess config file generate successfully at %s !" % (rawdata_dir + '/data_assess.cfg')
	
def make_dir(wanted_dir):
	if os.path.exists( wanted_dir ):
		print "%s is exists!" % wanted_dir
	#	os.system('rm -r ' + wanted_dir)
	#	os.mkdir(wanted_dir)
	else:
		os.mkdir(wanted_dir)	
	
def __main__():
	description = "This script is used to do data assess!\n"
	newParser = argparse.ArgumentParser( description = description );
	newParser.add_argument( "-rawdata", dest="rawdata");
	newParser.add_argument( "-mapping", dest="mapping");
	newParser.add_argument( "-ps", dest="pair_single_flag", help="paired/single");
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	rawdata_dir = argsDict['rawdata']
	mapping_file_url = argsDict['mapping']
	pair_single_flag = argsDict['pair_single_flag']
	
	start_time = datetime.datetime.now()
	
	make_dir(rawdata_dir + '/DataAssess')
	if pair_single_flag == 'paired':
		generate_data_assess_cfg(rawdata_dir, mapping_file_url)
		#perl /share/nas1/liuw/research/pipelines/Reseq_process_v1.5_lw/Tools/Data_Assess/saturation_random_stat_v1.1.pl -dataconfig data.config -Q 33 -o /share/nas1/liuw/project/BMK-Exom-Project/DataAssess
		subprocess.call('perl '+sys.path[0]+'/Data_Assess/saturation_random_stat_v1.1.pl -dataconfig ' + rawdata_dir + '/data_assess.cfg ' + ' -Q 33 -o ' + rawdata_dir + '/DataAssess', stdout = subprocess.PIPE, shell = True)
	else:
		generate_data_assess_cfg1(rawdata_dir, mapping_file_url)
		subprocess.call('perl '+sys.path[0]+'/Data_Assess/saturation_random_stat_v1.1_SE.pl -dataconfig ' + rawdata_dir + '/data_assess.cfg ' + ' -Q 33 -o ' + rawdata_dir + '/DataAssess', stdout = subprocess.PIPE, shell = True)
	end_time = datetime.datetime.now()
	print "Time Elapse:\t\t%s" % ( end_time-start_time ) 

if __name__ == "__main__": __main__()
