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
	newParser.add_argument( "-db", dest="db", help="blast database", default='/share/home/big/database/rpsH/rpsH');
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	seq = argsDict['seq']
	cpu = argsDict['cpu']
	db = argsDict['db']
	outdir = argsDict['outdir']
	
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	pick_otu_sh = outdir + '/pick_otus.sh'
	biom = outdir+'/picked_otus_uclust/otu_table_mc2_w_tax.biom'
	os.system('mkdir -p '+outdir+'/tmp')
	os.system('rm -rf '+outdir+'/picked_otus_uclust')
	pick_otu_cmd_list = ['pick_otus.py -i '+seq+' -o '+outdir+'/picked_otus_uclust']
	#pick_otu_cmd_list = ['pick_de_novo_otus.py -i '+seq+' -p '+sys.path[0]+'/ucrss_params_denovo_for_genes.txt -aO 6 -o '+outdir+'/picked_otus_uclust']
	pick_otu_cmd_list.append('pick_rep_set.py -i '+outdir +'/picked_otus_uclust/final_seqs_otus.txt -f '+seq+' -m most_abundant -o '+outdir+'/picked_otus_uclust/rep_set.fna')
	#pick_otu_cmd_list.append('blastx -query '+outdir+'/picked_otus_uclust/rep_set.fna -db '+db+' -evalue 1e-5 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 10 -outfmt 5 -num_threads '+str(cpu)+' -out '+outdir+'/picked_otus_uclust/rep_set_blast.result')
	pick_otu_cmd_list.append('blastx -query '+outdir+'/picked_otus_uclust/rep_set.fna -db '+db+' -evalue 1e-5 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 10 -outfmt 5 -out '+outdir+'/picked_otus_uclust/rep_set_blast.result')
	#get taxa assign result: taxa_assign.txt
	pick_otu_cmd_list.append('python '+sys.path[0]+'/biopython_xml_for_nr_blast.py '+outdir+'/picked_otus_uclust/rep_set_blast.result 0.7')
	pick_otu_cmd_list.append('make_otu_table.py -i '+outdir+'/picked_otus_uclust/final_seqs_otus.txt -t '+outdir+'/picked_otus_uclust/taxa_assign.txt -o '+biom)
	for cmd in pick_otu_cmd_list:
		os.system(cmd)
	os.system('biom summarize-table -i '+biom+' -o '+outdir+'/seqs_per_sample.txt')
	os.system('biom summarize-table -i '+biom+' --qualitative -o '+outdir+'/otus_per_sample.txt')
	os.system('mv '+biom+' '+outdir)
	os.system('mv '+outdir+'/picked_otus_uclust/rep_set.fna '+outdir)
	#os.system('rm -rf '+outdir+'/picked_otus_uclust')
	os.system('rm -rf '+outdir+'/tmp')
	
	pick_otu_handle = open(pick_otu_sh, 'w')
	cmd_list = []
	cmd_list.extend(pick_otu_cmd_list)
	cmd_list.append('biom summarize-table -i '+biom+' -o '+outdir+'/seqs_per_sample.txt')
	cmd_list.append('biom summarize-table -i '+biom+' --qualitative -o '+outdir+'/otus_per_sample.txt')
	cmd_list.append('mv '+biom+' '+outdir)
	cmd_list.append('mv '+outdir+'/picked_otus_uclust/rep_set.fna '+outdir)
	cmd_list.append('rm -rf '+outdir+'/picked_otus_uclust')	
	pick_otu_handle.write('\n'.join(cmd_list))
	pick_otu_handle.close()
		
if __name__ == "__main__": __main__()