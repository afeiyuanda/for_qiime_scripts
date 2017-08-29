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
-mapping mapping.txt 
-cpu 4
-outdir /path/to/outdir
		"""

def get_groups(mapping_file):
	f = open(mapping_file)
	lines = f.readlines()
	f.close()
	names = ['#SampleID','BarcodeSequence','LinkerPrimerSequence','ReversePrimer','Description']
	groups_list = []
	title_list = lines[0].strip().split('\t')
	for title in title_list:
		if title not in names:
			groups_list.append(title)
	return groups_list

def read_seqs_file(seqs_per_sample_file,start_line_flag):
	#seqs_per_sample_file:start_line_flag is 16
	#otus_per_sample_file:start_line_flag is 14
	f = open( seqs_per_sample_file )
	lines = f.readlines()
	f.close()
	d = {}
	for line in lines[int(start_line_flag):]:
		line_list = line.rstrip().split(' ')
		sampleid = line_list[1][:-1]
		seqs_num = line_list[2]
		d[sampleid] = seqs_num
	print d
	return d
	
def __main__():

	newParser = argparse.ArgumentParser( description = usage() );
	newParser.add_argument( "-seq", dest="seq", help="final_seqs.fna", required=True);
	newParser.add_argument( "-mapping", dest="mapping", help="mapping files containing group infor", required=True);
	#newParser.add_argument( "-biom", dest="biom", help="biom file", required=True);
	#newParser.add_argument( "-tre", dest="tre", help="tre file", required=True);
	newParser.add_argument( "-cpu", dest="cpu", help="cpu number", default='2');
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	seq = argsDict['seq']
	mapping = argsDict['mapping']
	#biom = argsDict['biom']
	#tre = argsDict['tre']
	cpu = argsDict['cpu']
	outdir = argsDict['outdir']
	
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	
	groups_list = get_groups(mapping)
	pick_otu_sh = outdir + '/pick_otu.sh'
	pick_otu_handle = open(pick_otu_sh, 'w')
	pick_cmd = 'pick_open_reference_otus.py -i '+seq+' -r $reference_seqs  -o '+outdir+'/picked_otus_uclust -p '+sys.path[0]+'/ucrss_params.txt -m uclust -aO '+str(cpu)+'\n'
	pick_otu_handle.write('rm -rf '+outdir+'/picked_otus_uclust\n')
	pick_otu_handle.write(pick_cmd)
	pick_otu_handle.write('biom  summarize-table -i '+outdir+'/picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom -o '+outdir+'/seqs_per_sample.txt\n' )
	pick_otu_handle.write('biom  summarize-table -i '+outdir+'/picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom --qualitative -o '+outdir+'/otus_per_sample.txt\n' )
	#subprocess.call( 'sh '+pick_otu_sh, stdout=subprocess.PIPE, shell=True )
	seqs_per_sample_dic = read_seqs_file(outdir+'/seqs_per_sample.txt', 16)
	seqs_num_list = [float(i) for i in seqs_per_sample_dic.values()]
	min_seqs_num = min(seqs_num_list)
	max_seqs_num = max(seqs_num_list)
	
	for group in groups_list:
		#Alpha diversity analyse
		alpha_sh = outdir + '/alpha_'+group+'.sh'
		alpha_handle = open(alpha_sh, 'w')
		alpha_handle.write('sort_otu_table.py -i '+outdir+'/picked_otus_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom -m '+mapping+' -s '+group+' -o '+outdir+'/sorted_otu_table_'+group+'.biom\n')
		m = int(min(1000, min_seqs_num))
		x = int(max_seqs_num - 100)
		s = int((x-m)/20)
		alpha_handle.write('parallel_multiple_rarefactions.py -i '+outdir+'/sorted_otu_table_'+group+'.biom -m '+str(m)+' -x '+str(x)+' -s '+str(s)+' -n 2 -o '+outdir+'/alpha_diversity_'+group+'/rarefied_otu_tables -O '+str(cpu)+'\n')
		alpha_handle.write('parallel_alpha_diversity.py -i '+outdir+'/alpha_diversity_'+group+'/rarefied_otu_tables -m ACE,simpson,shannon,PD_whole_tree,chao1,observed_species,goods_coverage -o '+outdir+'/alpha_diversity_'+group+'/adiv/ -t '+outdir+'/picked_otus_uclust/rep_set.tre -O '+str(cpu)+'\n')
		alpha_handle.write('collate_alpha.py -i '+outdir+'/alpha_diversity_'+group+'/adiv/ -o '+outdir+'/alpha_diversity_'+group+'/collated_alpha/\n')
		alpha_handle.write('python '+sys.path[0]+'/plot_rarefaction_curve.py '+outdir+'/alpha_diversity_'+group+'/collated_alpha/ bottomright\n')
		alpha_handle.write('python '+sys.path[0]+'/boxplot_rarefaction_curve.py '+outdir+'/alpha_diversity_'+group+'/collated_alpha/ '+mapping+' '+group+' bottomright\n')
		alpha_handle.close()

		#Beta diversity analyse
		beta_sh = outdir + '/beta_'+group+'.sh'
		beta_handle = open(beta_sh, 'w')
		beta_handle.write('mkdir -p '+outdir+'/beta_diversity\n')
		beta_handle.write('single_rarefaction.py -i '+outdir+'/sorted_otu_table_'+group+'.biom -o '+outdir+'/beta_diversity/otu_table_even.biom -d '+str(int(min_seqs_num*0.95))+'\n')
		beta_handle.write('parallel_beta_diversity.py -i '+outdir+'/beta_diversity/otu_table_even.biom -t '+outdir+'/picked_otus_uclust/rep_set.tre -o '+outdir+'/beta_diversity/beta_div -O '+str(4)+'\n')
		beta_handle.write('principal_coordinates.py -i '+outdir+'/beta_diversity/beta_div/ -o '+outdir+'/beta_diversity/beta_div_pca_result\n')
		beta_handle.write('python '+sys.path[0]+'/plot_pca_with_pc_vectors.py '+outdir+'/beta_diversity/beta_div_pca_result/pcoa_weighted_unifrac_otu_table_even.txt '+mapping+' '+group+' bottomright\n')
		beta_handle.write('python '+sys.path[0]+'/plot_pca_with_pc_vectors.py '+outdir+'/beta_diversity/beta_div_pca_result/pcoa_unweighted_unifrac_otu_table_even.txt '+mapping+' '+group+' bottomright\n')
		beta_handle.write('biom  summarize-table -i '+outdir+'/beta_diversity/otu_table_even.biom -o '+outdir+'/even_seqs_per_sample.txt\n' )
		beta_handle.write('biom  summarize-table -i '+outdir+'/beta_diversity/otu_table_even.biom --qualitative -o '+outdir+'/even_otus_per_sample.txt\n' )
		beta_handle.write('alpha_diversity.py -i '+outdir+'/beta_diversity/otu_table_even.biom -m ACE,simpson,shannon,PD_whole_tree,chao1,observed_species,goods_coverage -o '+outdir+'/even_adiv.txt -t '+outdir+'/picked_otus_uclust/rep_set.tre\n')
		beta_handle.write('summarize_otu_by_cat.py -i '+outdir+'/beta_diversity/otu_table_even.biom   -m '+mapping+' -c '+group+' -o '+outdir+'/even_sorted_otu_table_by_'+group+'.biom\n')
		beta_handle.write('plot_rank_abundance_graph.py -i '+outdir+'/even_sorted_otu_table_by_'+group+'.biom -s \'*\' -x -v -f png -o '+outdir+'/OTU_rank_abundance_even.png\n')
		beta_handle.write('summarize_taxa_through_plots.py -i '+outdir+'/beta_diversity/otu_table_even.biom -o '+outdir+'/wf_taxa_summary_even -m '+mapping+'\n')
		beta_handle.write('jackknifed_beta_diversity.py -i '+outdir+'/sorted_otu_table_'+group+'.biom -t '+outdir+'/picked_otus_uclust/rep_set.tre  -m '+mapping+' -o '+outdir+'/wf_jack -e '+str(int(min_seqs_num*0.95))+'\n')
		beta_handle.write('python '+sys.path[0]+'/generate_upgma_tree_images.py '+outdir+'\n')
		beta_handle.close()
		
		#Ploting
		plot_sh = outdir + '/plot_'+group+'.sh'
		plot_handle = open(plot_sh, 'w')
		plot_handle.write('mkdir -p '+outdir+'/OTU_Heatmap_even\n')
		plot_handle.write('biom convert -i '+outdir+'/beta_diversity/otu_table_even.biom  -b  -o '+outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name_'+group+'.txt --header-key taxonomy\n')
		plot_handle.write('python '+sys.path[0]+'/plot_heatmap.py  '+outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name_'+group+'.txt 80 TRUE\n')
		plot_handle.write('mkdir -p '+outdir+'/OTU_Venn_diagram_even\n')
		plot_handle.write('biom convert -i '+outdir+'/even_sorted_otu_table_by_'+group+'.biom  -b  -o '+outdir+'/OTU_Venn_diagram_even/even_sorted_otu_table_by_'+group+'.txt\n')
		plot_handle.write('python '+sys.path[0]+'/plot_venn_diagram.py  '+outdir+'/OTU_Venn_diagram_even/even_sorted_otu_table_by_'+group+'.txt\n')
		
		plot_handle.write('python '+sys.path[0]+'/taxo_barplot_for30.56.py '+outdir+'/wf_taxa_summary_even/  '+outdir+'/taxa_summary_even\n')
		plot_handle.write('python '+sys.path[0]+'/plot_average_taxo_abundance.py '+outdir+'/taxa_summary_even/ '+mapping+' '+group+'\n')
		plot_handle.write('python '+sys.path[0]+'/taxonomy_anova.py '+outdir+'/taxa_summary_even/ '+mapping+' '+group+'\n')
		plot_handle.write('python '+sys.path[0]+'/collated_seqs_otus_num.py '+mapping+' '+group+' '+outdir+'/seqs_per_sample.txt '+outdir+'/otus_per_sample.txt '+outdir+'/even_seqs_per_sample.txt '+outdir+'/even_otus_per_sample.txt '+outdir+'/total_summary.txt\n')

		plot_handle.close()
		final_result_dir = outdir+'/Analyse_Resluts_'+group
		os.system('mkdir '+final_result_dir)
		os.system('mkdir -p '+final_result_dir+'/Alpha_diversity')
		os.system('mkdir -p '+final_result_dir+'/OTU_Heatmap_even')
		os.system('mkdir -p '+final_result_dir+'/OTU_Venn_diagram_even')
		os.system('mkdir -p '+final_result_dir+'/PCoA')
		os.system('mkdir -p '+final_result_dir+'/PCoA/3DPCoA_weighted_unifrac')
		os.system('mkdir -p '+final_result_dir+'/PCoA/3DPCoA_unweighted_unifrac')
		os.system('mkdir -p '+final_result_dir+'/UPGMA')
		os.system('mkdir -p '+final_result_dir+'/taxa_summary_even')
		
		os.system('mv '+outdir+'/total_summary.txt '+final_result_dir)
		os.system('mv '+outdir+'/even_adiv.txt '+final_result_dir)
		os.system('mv '+outdir+'/OTU_Heatmap_even '+final_result_dir)
		os.system('rm '+final_result_dir+'/OTU_Heatmap_even/*.r*')
		os.system('mv '+outdir+'/OTU_Venn_diagram_even '+final_result_dir)
		os.system('rm '+final_result_dir+'/OTU_Venn_diagram_even/*.r*')
		os.system('mv '+outdir+'/alpha_diversity_'+group+'/collated_alpha/*.txt '+final_result_dir+'/Alpha_diversity')
		os.system('mv '+outdir+'/alpha_diversity_'+group+'/collated_alpha/*.png '+final_result_dir+'/Alpha_diversity')
		os.system('mv '+outdir+'/beta_diversity/beta_div_pca_result/*.png '+final_result_dir+'/PCoA')
		os.system('mv '+outdir+'/wf_jack/weighted_unifrac/emperor_pcoa_plots/* '+final_result_dir+'/PCoA/3DPCoA_weighted_unifrac')
		os.system('mv '+outdir+'/wf_jack/unweighted_unifrac/emperor_pcoa_plots/* '+final_result_dir+'/PCoA/3DPCoA_unweighted_unifrac')
		os.system('mv '+outdir+'/taxa_summary_even '+final_result_dir)
		os.system('rm '+final_result_dir+'/taxa_summary_even/*/*.r*')
		os.system('mv '+outdir+'/UPGMA '+final_result_dir)
		os.system('mv '+outdir+'/OTU_rank_abundance_even.png '+final_result_dir)
		
if __name__ == "__main__": __main__()