#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.15
'''

import sys,os,glob,subprocess,re,argparse

def usage():
    print """
Usage:	python """+sys.argv[0]+ """ 
-seq seq 
-cpu 4
-outdir /path/to/outdir
		"""
reference_seqs = '/share/nas2/genome/biosoft/QIIME/gg_otus-13_8-release/rep_set/97_otus.fasta'
	
def __main__():

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-seq", dest="seq", help="final_seqs.fna", required=True);
	newParser.add_argument( "-cpu", dest="cpu", help="cpu number", default='2');
	#newParser.add_argument( "-queue", dest="queue", help="queue name", default='bioloong');
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	seq = argsDict['seq']
	cpu = argsDict['cpu']
	#queue = argsDict['queue']
	outdir = argsDict['outdir']
	
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	pick_otu_sh = outdir + '/pick_otus.sh'
	biom = outdir+'/picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom'
	
	os.system('rm -rf '+outdir+'/picked_otus_uclust')
	os.system('pick_open_reference_otus.py -i '+seq+' -r /share/nas2/genome/biosoft/QIIME/gg_otus-13_8-release/rep_set/97_otus.fasta  -o '+outdir+'/picked_otus_uclust -p '+sys.path[0]+'/ucrss_params.txt -m uclust -aO '+str(cpu))
	os.system('biom summarize-table -i '+biom+' -o '+outdir+'/seqs_per_sample.txt')
	os.system('biom summarize-table -i '+biom+' --qualitative -o '+outdir+'/otus_per_sample.txt')
	os.system('mv '+outdir+'/picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom '+outdir)
	os.system('mv '+outdir+'/picked_otus_uclust/rep_set.tre '+outdir)
	os.system('mv '+outdir+'/picked_otus_uclust/rep_set.fna '+outdir)
	#os.system('rm -rf '+outdir+'/picked_otus_uclust')
	os.system('rm -rf '+outdir+'/tmp')
	
	pick_otu_handle = open(pick_otu_sh, 'w')
	cmd_list = []
	cmd_list.append('pick_open_reference_otus.py -i '+seq+' -r /share/nas2/genome/biosoft/QIIME/gg_otus-13_8-release/rep_set/97_otus.fasta  -o '+outdir+'/picked_otus_uclust -p '+sys.path[0]+'/ucrss_params.txt -m uclust -aO '+str(cpu))
	cmd_list.append('biom summarize-table -i '+biom+' -o '+outdir+'/seqs_per_sample.txt')
	cmd_list.append('biom summarize-table -i '+biom+' --qualitative -o '+outdir+'/otus_per_sample.txt')
	cmd_list.append('mv '+outdir+'/picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom '+outdir)
	cmd_list.append('mv '+outdir+'/picked_otus_uclust/rep_set.tre '+outdir)
	cmd_list.append('mv '+outdir+'/picked_otus_uclust/rep_set.fna '+outdir)
	cmd_list.append('rm -rf '+outdir+'/picked_otus_uclust')	
	pick_otu_handle.write('\n'.join(cmd_list))
	pick_otu_handle.close()
		
if __name__ == "__main__": __main__()
