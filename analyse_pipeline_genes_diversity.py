#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.12.15
'''

import sys,os,glob,subprocess,re,argparse,string

def usage():
    print """
Usage:	python """+sys.argv[0]+ """ 
-mapping /share/home/big/liuwei_test/mapping.txt 
-biom /share/home/big/liuwei_test/picked_otus_uclust/otu_table_mc2_w_tax.biom 
-group Treatment 
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
	newParser.add_argument( "-mapping", dest="mapping", help="mapping files containing group infor", required=True);
	newParser.add_argument( "-biom", dest="biom", help="biom file", required=True);
	newParser.add_argument( "-group", dest="gp", help="group name in mapping file");
	newParser.add_argument( "-cpu", dest="cpu", help="cpu number", default='2');
	newParser.add_argument( "-norm", dest="norm", help="normalization");
	newParser.add_argument( "-test", dest="test", help="significant test:t_test,annova,wilcox", default='t_test');
	newParser.add_argument( "-outdir", dest="outdir", help="Your analysis output dir", default=os.getcwd() );
	
	
	args = newParser.parse_args();
	argsDict = args.__dict__;
	
	mapping = argsDict['mapping']
	biom = argsDict['biom']
	gp = argsDict['gp']
	cpu = argsDict['cpu']
	norm = argsDict['norm']
	test = argsDict['test']
	outdir = argsDict['outdir']
	
	if outdir != os.getcwd():
		if os.path.exists( outdir ):
			pass
		else:
			os.makedirs( outdir )
	biom_file_dir = os.path.dirname(biom)
	#print biom_file_dir
	if os.path.exists(biom_file_dir+'/seqs_per_sample.txt'):
		os.system('rm '+biom_file_dir+'/seqs_per_sample.txt')
		os.system('rm '+biom_file_dir+'/otus_per_sample.txt')
		os.system('biom summarize-table -i '+biom+' -o '+biom_file_dir+'/seqs_per_sample.txt')
		os.system('biom summarize-table -i '+biom+' --qualitative -o '+biom_file_dir+'/otus_per_sample.txt')
		seqs_per_sample_dic = read_seqs_file(biom_file_dir+'/seqs_per_sample.txt', 16)
	else:
		os.system('biom summarize-table -i '+biom+' -o '+biom_file_dir+'/seqs_per_sample.txt')
		os.system('biom summarize-table -i '+biom+' --qualitative -o '+biom_file_dir+'/otus_per_sample.txt')
		seqs_per_sample_dic = read_seqs_file(biom_file_dir+'/seqs_per_sample.txt', 16)
	seqs_num_list = sorted([float(i) for i in seqs_per_sample_dic.values()])
	reverse_seqs_num_list = seqs_num_list[::-1]
	min_seqs_num = min(seqs_num_list)
	#max_seqs_num = max(seqs_num_list)
	if min_seqs_num == 1:
		min_seqs_num = 100
	m = max(500, int(min_seqs_num*0.2))
	for i in range(len(seqs_num_list)):
		if (reverse_seqs_num_list[i]-reverse_seqs_num_list[i+1])/reverse_seqs_num_list[i+1] > 0.7:
			continue
		else:
			max_seqs_num = int(reverse_seqs_num_list[i])
			break
	print 'max_seqs_num is %s' % max_seqs_num
	x = max_seqs_num
	s = int((max_seqs_num-m)/20)	
	
	if norm is None:
		norm = min_seqs_num
	depth = int(min_seqs_num*0.95)
	
	if gp is None:
		group_list = get_groups(mapping)
	else:
		group_list = [gp]
	
	for group in group_list:
		new_outdir = outdir + '/'+string.lower(group)+'_analysis'
		if os.path.exists( new_outdir ):
			#os.system('rm -rf new_outdir')
			pass
		else:
			os.makedirs( new_outdir )
		all_cmd_list = []
		all_cmd_list.append('sort_otu_table.py -i '+biom+' -m '+mapping+' -s '+group+' -o '+new_outdir+'/sorted_otu_table.biom')
		all_cmd_list.append('parallel_multiple_rarefactions.py -i '+new_outdir+'/sorted_otu_table.biom -m '+str(m)+' -x '+str(max_seqs_num)+' -s '+str(s)+' -n '+str(cpu)+' -o '+new_outdir+'/alpha_diversity/rarefied_otu_tables -O '+str(cpu))
		#all_cmd_list.append('parallel_alpha_diversity.py -i '+new_outdir+'/alpha_diversity/rarefied_otu_tables -m ACE,simpson,shannon,chao1,observed_species,goods_coverage -o '+new_outdir+'/alpha_diversity/adiv/ -O '+str(cpu))
		all_cmd_list.append('parallel_alpha_diversity.py -i '+new_outdir+'/alpha_diversity/rarefied_otu_tables -m simpson,shannon,chao1,observed_species,goods_coverage -o '+new_outdir+'/alpha_diversity/adiv/ -O '+str(cpu))
		all_cmd_list.append('collate_alpha.py -i '+new_outdir+'/alpha_diversity/adiv/ -o '+new_outdir+'/alpha_diversity/collated_alpha/')
		all_cmd_list.append('python '+sys.path[0]+'/plot_rarefaction_curve.py '+new_outdir+'/alpha_diversity/collated_alpha/ bottomright')
		all_cmd_list.append('/share/home/big/qiime_1.8/python-2.7.3-release/bin/python '+sys.path[0]+'/plot_rarefaction_curve_with_error_bar.py '+new_outdir+'/alpha_diversity/collated_alpha/ '+mapping+' '+group)
		#all_cmd_list.append('python '+sys.path[0]+'/boxplot_rarefaction_curve.py '+new_outdir+'/alpha_diversity/collated_alpha/ '+mapping+' '+group+' bottomright')
		all_cmd_list.append('mkdir -p '+new_outdir+'/beta_diversity')
		all_cmd_list.append('single_rarefaction.py -i '+new_outdir+'/sorted_otu_table.biom -o '+new_outdir+'/beta_diversity/otu_table_even.biom -d '+str(depth))
		all_cmd_list.append('multiple_rarefactions_even_depth.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom -d '+str(int(depth*0.9))+' -o '+new_outdir+'/beta_diversity/rarefaction/')
		all_cmd_list.append('beta_diversity.py -i '+new_outdir+'/beta_diversity/rarefaction/ -o '+new_outdir+'/beta_diversity/rare_dm/ -m bray_curtis')
		all_cmd_list.append('upgma_cluster.py -i '+new_outdir+'/beta_diversity/rare_dm/ -o '+new_outdir+'/beta_diversity/rare_upgma/')
		all_cmd_list.append('consensus_tree.py -i '+new_outdir+'/beta_diversity/rare_upgma -o '+new_outdir+'/beta_diversity/rare_upgma_consensus.tre')
		all_cmd_list.append('tree_compare.py -s '+new_outdir+'/beta_diversity/rare_upgma -m '+new_outdir+'/beta_diversity/rare_upgma_consensus.tre -o '+new_outdir+'/beta_diversity/upgma_cmp/')
		
		all_cmd_list.append('biom  summarize-table -i '+new_outdir+'/beta_diversity/otu_table_even.biom -o '+new_outdir+'/even_seqs_per_sample.txt')
		all_cmd_list.append('biom  summarize-table -i '+new_outdir+'/beta_diversity/otu_table_even.biom --qualitative -o '+new_outdir+'/even_otus_per_sample.txt')
		all_cmd_list.append('alpha_diversity.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom -m ACE,simpson,shannon,chao1,observed_species,goods_coverage -o '+new_outdir+'/even_adiv.txt')
		#OTU analyse
		all_cmd_list.append('summarize_otu_by_cat.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom   -m '+mapping+' -c '+group+' -o '+new_outdir+'/even_sorted_otu_table_by_'+group+'.biom')
		all_cmd_list.append('plot_rank_abundance_graph.py -i '+new_outdir+'/even_sorted_otu_table_by_'+group+'.biom -s \'*\' -x -v -f png -o '+new_outdir+'/OTU_rank_abundance_even.png')
		all_cmd_list.append('mkdir -p '+new_outdir+'/OTU_Heatmap_even')
		all_cmd_list.append('biom convert -i '+new_outdir+'/beta_diversity/otu_table_even.biom  -b  -o '+new_outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name.txt --header-key taxonomy')
		all_cmd_list.append('python '+sys.path[0]+'/plot_heatmap_enhanced.py  '+new_outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name.txt 80 TRUE otu')
		all_cmd_list.append('python '+sys.path[0]+'/plot_heatmap_enhanced.py  '+new_outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name.txt 80 TRUE taxonomy')
		all_cmd_list.append('mkdir -p '+new_outdir+'/OTU_Venn_diagram_even')
		all_cmd_list.append('biom convert -i '+new_outdir+'/even_sorted_otu_table_by_'+group+'.biom  -b  -o '+new_outdir+'/OTU_Venn_diagram_even/even_sorted_otu_table_by_'+group+'.txt')
		all_cmd_list.append('python '+sys.path[0]+'/plot_venn_diagram.py  '+new_outdir+'/OTU_Venn_diagram_even/even_sorted_otu_table_by_'+group+'.txt')
		#summarize_taxa
		all_cmd_list.append('summarize_taxa.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom -L 1 -o '+new_outdir+'/wf_taxa_summary_even')		
		all_cmd_list.append('python '+sys.path[0]+'/taxo_barplot_for_only_species.py '+new_outdir+'/wf_taxa_summary_even/  '+new_outdir+'/taxa_summary_even')
		all_cmd_list.append('python '+sys.path[0]+'/plot_average_taxo_abundance.py '+new_outdir+'/taxa_summary_even/ '+mapping+' '+group)
		#beta_diversity analyse
		all_cmd_list.append('beta_diversity_through_plots.py -i '+biom+' -o '+new_outdir+'/bdiv_even_bray_curtis -e '+str(depth)+' -m '+mapping+' -p '+sys.path[0]+'/its_beta_diversity_params.txt')
		all_cmd_list.append('make_3d_plots.py -i '+new_outdir+'/bdiv_even_bray_curtis/bray_curtis_pc.txt '+' -m '+mapping+' -o '+new_outdir+'/plot3d_bray_curtis')
		all_cmd_list.append('python '+sys.path[0]+'/plot_pca_with_pc_vectors.py '+new_outdir+'/bdiv_even_bray_curtis/bray_curtis_pc.txt '+mapping+' '+group+' bottomright')

		os.system('sort_otu_table.py -i '+biom+' -m '+mapping+' -s '+group+' -o '+new_outdir+'/sorted_otu_table.biom')
		os.system('parallel_multiple_rarefactions.py -i '+new_outdir+'/sorted_otu_table.biom -m '+str(m)+' -x '+str(max_seqs_num)+' -s '+str(s)+' -n '+str(cpu)+' -o '+new_outdir+'/alpha_diversity/rarefied_otu_tables -O '+str(cpu))
		os.system('parallel_alpha_diversity.py -i '+new_outdir+'/alpha_diversity/rarefied_otu_tables -m simpson,shannon,chao1,observed_species,goods_coverage -o '+new_outdir+'/alpha_diversity/adiv/ -O '+str(cpu))
		#os.system('parallel_alpha_diversity.py -i '+new_outdir+'/alpha_diversity/rarefied_otu_tables -m ACE,simpson,shannon,chao1,observed_species,goods_coverage -o '+new_outdir+'/alpha_diversity/adiv/ -O '+str(cpu))
		os.system('collate_alpha.py -i '+new_outdir+'/alpha_diversity/adiv/ -o '+new_outdir+'/alpha_diversity/collated_alpha/')
		os.system('python '+sys.path[0]+'/plot_rarefaction_curve.py '+new_outdir+'/alpha_diversity/collated_alpha/ bottomright')
		os.system('/share/home/big/qiime_1.8/python-2.7.3-release/bin/python '+sys.path[0]+'/plot_rarefaction_curve_with_error_bar.py '+new_outdir+'/alpha_diversity/collated_alpha/ '+mapping+' '+group)
		#os.system('python '+sys.path[0]+'/boxplot_rarefaction_curve.py '+new_outdir+'/alpha_diversity/collated_alpha/ '+mapping+' '+group+' bottomright')
		os.system('mkdir -p '+new_outdir+'/beta_diversity')
		os.system('single_rarefaction.py -i '+new_outdir+'/sorted_otu_table.biom -o '+new_outdir+'/beta_diversity/otu_table_even.biom -d '+str(depth))
		os.system('multiple_rarefactions_even_depth.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom -d '+str(int(depth*0.9))+' -o '+new_outdir+'/beta_diversity/rarefaction/')
		os.system('beta_diversity.py -i '+new_outdir+'/beta_diversity/rarefaction/ -o '+new_outdir+'/beta_diversity/rare_dm/ -m bray_curtis')
		os.system('upgma_cluster.py -i '+new_outdir+'/beta_diversity/rare_dm/ -o '+new_outdir+'/beta_diversity/rare_upgma/')
		os.system('consensus_tree.py -i '+new_outdir+'/beta_diversity/rare_upgma -o '+new_outdir+'/beta_diversity/rare_upgma_consensus.tre')
		os.system('tree_compare.py -s '+new_outdir+'/beta_diversity/rare_upgma -m '+new_outdir+'/beta_diversity/rare_upgma_consensus.tre -o '+new_outdir+'/beta_diversity/upgma_cmp/')
		
		os.system('biom  summarize-table -i '+new_outdir+'/beta_diversity/otu_table_even.biom -o '+new_outdir+'/even_seqs_per_sample.txt')
		os.system('biom  summarize-table -i '+new_outdir+'/beta_diversity/otu_table_even.biom --qualitative -o '+new_outdir+'/even_otus_per_sample.txt')
		os.system('alpha_diversity.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom -m ACE,simpson,shannon,chao1,observed_species,goods_coverage -o '+new_outdir+'/even_adiv.txt')
		#OTU analyse
		os.system('summarize_otu_by_cat.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom   -m '+mapping+' -c '+group+' -o '+new_outdir+'/even_sorted_otu_table_by_'+group+'.biom')
		os.system('plot_rank_abundance_graph.py -i '+new_outdir+'/even_sorted_otu_table_by_'+group+'.biom -s \'*\' -x -v -f png -o '+new_outdir+'/OTU_rank_abundance_even.png')
		os.system('mkdir -p '+new_outdir+'/OTU_Heatmap_even')
		os.system('biom convert -i '+new_outdir+'/beta_diversity/otu_table_even.biom  -b  -o '+new_outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name.txt --header-key taxonomy')
		os.system('python '+sys.path[0]+'/plot_heatmap_enhanced.py  '+new_outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name.txt 80 TRUE otu')
		os.system('python '+sys.path[0]+'/plot_heatmap_enhanced.py  '+new_outdir+'/OTU_Heatmap_even/even_otu_table_taxo_name.txt 80 TRUE taxonomy')
		os.system('mkdir -p '+new_outdir+'/OTU_Venn_diagram_even')
		os.system('biom convert -i '+new_outdir+'/even_sorted_otu_table_by_'+group+'.biom  -b  -o '+new_outdir+'/OTU_Venn_diagram_even/even_sorted_otu_table_by_'+group+'.txt')
		os.system('python '+sys.path[0]+'/plot_venn_diagram.py  '+new_outdir+'/OTU_Venn_diagram_even/even_sorted_otu_table_by_'+group+'.txt')
		#summarize_taxa
		os.system('summarize_taxa.py -i '+new_outdir+'/beta_diversity/otu_table_even.biom -L 1 -o '+new_outdir+'/wf_taxa_summary_even')		
		os.system('python '+sys.path[0]+'/taxo_barplot_for_only_species.py '+new_outdir+'/wf_taxa_summary_even/  '+new_outdir+'/taxa_summary_even')
		os.system('python '+sys.path[0]+'/plot_average_taxo_abundance.py '+new_outdir+'/taxa_summary_even/ '+mapping+' '+group)
		#beta_diversity analyse
		os.system('beta_diversity_through_plots.py -i '+biom+' -o '+new_outdir+'/bdiv_even_bray_curtis -e '+str(depth)+' -m '+mapping+' -p '+sys.path[0]+'/its_beta_diversity_params.txt')
		os.system('make_3d_plots.py -i '+new_outdir+'/bdiv_even_bray_curtis/bray_curtis_pc.txt '+' -m '+mapping+' -o '+new_outdir+'/plot3d_bray_curtis')
		os.system('python '+sys.path[0]+'/plot_pca_with_pc_vectors.py '+new_outdir+'/bdiv_even_bray_curtis/bray_curtis_pc.txt '+mapping+' '+group+' bottomright')
	
		
		if test == 't_test':
			significant_cmd = 'python '+sys.path[0]+'/t_test.py '+new_outdir+'/taxa_summary_even/ '+mapping+' '+group
		elif test == 'annova':
			significant_cmd = 'python '+sys.path[0]+'/taxonomy_anova.py '+new_outdir+'/taxa_summary_even/ '+mapping+' '+group
		else:
			significant_cmd = 'python '+sys.path[0]+'/wilcox_test.py '+new_outdir+'/taxa_summary_even/ '+mapping+' '+group
		all_cmd_list.append(significant_cmd)
		all_cmd_list.append('python '+sys.path[0]+'/collated_seqs_otus_num.py '+mapping+' '+group+' '+biom_file_dir+'/seqs_per_sample.txt '+biom_file_dir+'/otus_per_sample.txt '+new_outdir+'/even_seqs_per_sample.txt '+new_outdir+'/even_otus_per_sample.txt '+new_outdir+'/total_summary.txt')
		os.system(significant_cmd)
		os.system('python '+sys.path[0]+'/collated_seqs_otus_num.py '+mapping+' '+group+' '+biom_file_dir+'/seqs_per_sample.txt '+biom_file_dir+'/otus_per_sample.txt '+new_outdir+'/even_seqs_per_sample.txt '+new_outdir+'/even_otus_per_sample.txt '+new_outdir+'/total_summary.txt')
		
		
		collate_results_cmd = []
		final_result_dir = outdir+'/'+string.capitalize(group)+'_Analyse_Resluts'
		collate_results_cmd.append('mkdir '+final_result_dir)
		collate_results_cmd.append('mkdir -p '+final_result_dir+'/Alpha_diversity')
		collate_results_cmd.append('mkdir -p '+final_result_dir+'/OTU_Heatmap_even')
		collate_results_cmd.append('mkdir -p '+final_result_dir+'/OTU_Venn_diagram_even')
		collate_results_cmd.append('mkdir -p '+final_result_dir+'/PCoA')
		collate_results_cmd.append('mkdir -p '+final_result_dir+'/taxa_summary_even')
		collate_results_cmd.append('mkdir -p '+final_result_dir+'/bray_curtis_distance')	
		collate_results_cmd.append('mv '+new_outdir+'/total_summary.txt '+final_result_dir)
		collate_results_cmd.append('mv '+new_outdir+'/even_adiv.txt '+final_result_dir)
		collate_results_cmd.append('mv '+new_outdir+'/OTU_Heatmap_even '+final_result_dir)
		collate_results_cmd.append('rm '+final_result_dir+'/OTU_Heatmap_even/*.r*')
		collate_results_cmd.append('mv '+new_outdir+'/OTU_Venn_diagram_even '+final_result_dir)
		collate_results_cmd.append('rm '+final_result_dir+'/OTU_Venn_diagram_even/*.r*')
		collate_results_cmd.append('mv '+new_outdir+'/alpha_diversity/collated_alpha/*.txt '+final_result_dir+'/Alpha_diversity')
		collate_results_cmd.append('mv '+new_outdir+'/alpha_diversity/collated_alpha/*.png '+final_result_dir+'/Alpha_diversity')
		collate_results_cmd.append('mv '+new_outdir+'/bdiv_even_bray_curtis/*.png '+final_result_dir+'/PCoA')
		collate_results_cmd.append('mv '+new_outdir+'/plot3d_bray_curtis '+final_result_dir+'/PCoA/3DPCoA_bray_curtis')
		collate_results_cmd.append('mv '+new_outdir+'/taxa_summary_even '+final_result_dir)
		collate_results_cmd.append('mv '+new_outdir+'/beta_diversity/upgma_cmp/ '+final_result_dir+'/UPGMA')
		collate_results_cmd.append('rm '+final_result_dir+'/taxa_summary_even/*/*.r*')
		collate_results_cmd.append('mv '+new_outdir+'/OTU_rank_abundance_even.png '+final_result_dir)
		collate_results_cmd.append('cp -r '+new_outdir+'/bdiv_even_bray_curtis/bray_curtis_dm.txt '+final_result_dir+'/bray_curtis_distance')

		os.system('mkdir '+final_result_dir)
		os.system('mkdir -p '+final_result_dir+'/Alpha_diversity')
		os.system('mkdir -p '+final_result_dir+'/OTU_Heatmap_even')
		os.system('mkdir -p '+final_result_dir+'/OTU_Venn_diagram_even')
		os.system('mkdir -p '+final_result_dir+'/PCoA')
		os.system('mkdir -p '+final_result_dir+'/taxa_summary_even')
		os.system('mkdir -p '+final_result_dir+'/bray_curtis_distance')	
		os.system('mv '+new_outdir+'/total_summary.txt '+final_result_dir)
		os.system('mv '+new_outdir+'/even_adiv.txt '+final_result_dir)
		os.system('mv '+new_outdir+'/OTU_Heatmap_even '+final_result_dir)
		os.system('rm '+final_result_dir+'/OTU_Heatmap_even/*.r*')
		os.system('mv '+new_outdir+'/OTU_Venn_diagram_even '+final_result_dir)
		os.system('rm '+final_result_dir+'/OTU_Venn_diagram_even/*.r*')
		os.system('mv '+new_outdir+'/alpha_diversity/collated_alpha/*.txt '+final_result_dir+'/Alpha_diversity')
		os.system('mv '+new_outdir+'/alpha_diversity/collated_alpha/*.png '+final_result_dir+'/Alpha_diversity')
		os.system('mv '+new_outdir+'/bdiv_even_bray_curtis/*.png '+final_result_dir+'/PCoA')
		os.system('mv '+new_outdir+'/plot3d_bray_curtis '+final_result_dir+'/PCoA/3DPCoA_bray_curtis')
		os.system('mv '+new_outdir+'/taxa_summary_even '+final_result_dir)
		os.system('mv '+new_outdir+'/beta_diversity/upgma_cmp/ '+final_result_dir+'/UPGMA')
		os.system('rm '+final_result_dir+'/taxa_summary_even/*/*.r*')
		os.system('mv '+new_outdir+'/OTU_rank_abundance_even.png '+final_result_dir)
		os.system('cp -r '+new_outdir+'/bdiv_even_bray_curtis/bray_curtis_dm.txt '+final_result_dir+'/bray_curtis_distance')
	
		
		f = open(new_outdir+'/diversity.sh', 'w')
		all_cmd_list.extend(collate_results_cmd)
		f.writelines('\n'.join(all_cmd_list)+'\n')
		f.close()
		
		os.system('rm '+outdir+'/*.r*')
		#os.system('rm '+outdir+'/*.txt')
		os.system('rm '+outdir+'/*.RData')
		os.system('rm '+outdir+'/*.Rout')
		os.system('rm -rf '+outdir+'/tmp')
		print '%s analyse finished!' % group
	
if __name__ == "__main__": __main__()