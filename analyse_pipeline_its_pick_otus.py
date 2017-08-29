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
	
def __main__():

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-seq", dest="seq", help="final_seqs.fna", required=True);
	newParser.add_argument( "-cpu", dest="cpu", help="cpu number", default='2');
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	seq = argsDict['seq']
	cpu = argsDict['cpu']
	outdir = argsDict['outdir']
	
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	pick_otu_sh = outdir + '/pick_otus.sh'
	biom = outdir+'/picked_otus_uclust/otu_table_mc2_w_tax.biom'
	
	os.system('rm -rf '+outdir+'/picked_otus_uclust')
	pick_otu_cmd = 'pick_open_reference_otus.py -i '+seq+' -r $its_reference_seqs -f -o '+outdir+'/picked_otus_uclust -p '+sys.path[0]+'/ucrss_params_its.txt --suppress_align_and_tree -m uclust -aO '+str(cpu)
	os.system(pick_otu_cmd)
	os.system('biom summarize-table -i '+biom+' -o '+outdir+'/seqs_per_sample.txt')
	os.system('biom summarize-table -i '+biom+' --qualitative -o '+outdir+'/otus_per_sample.txt')
	os.system('mv '+biom+' '+outdir)
	os.system('mv '+outdir+'/picked_otus_uclust/rep_set.fna '+outdir)
	#os.system('rm -rf '+outdir+'/picked_otus_uclust')
	os.system('rm -rf '+outdir+'/tmp')
	
	pick_otu_handle = open(pick_otu_sh, 'w')
	cmd_list = []
	cmd_list.append(pick_otu_cmd)
	cmd_list.append('biom summarize-table -i '+biom+' -o '+outdir+'/seqs_per_sample.txt')
	cmd_list.append('biom summarize-table -i '+biom+' --qualitative -o '+outdir+'/otus_per_sample.txt')
	cmd_list.append('mv '+biom+' '+outdir)
	cmd_list.append('mv '+outdir+'/picked_otus_uclust/rep_set.fna '+outdir)
	cmd_list.append('rm -rf '+outdir+'/picked_otus_uclust')	
	pick_otu_handle.write('\n'.join(cmd_list))
	pick_otu_handle.close()
		
if __name__ == "__main__": __main__()