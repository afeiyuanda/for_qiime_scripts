beta_diversity.py -i sorted_otu_table.biom -o wf_jack  -t picked_otus_uclust/rep_set.tre
multiple_rarefactions_even_depth.py -i sorted_otu_table.biom -d 9500 -o wf_jack/rarefaction/ 

mkdir wf_jack/weighted_unifrac/
upgma_cluster.py -i wf_jack/weighted_unifrac_sorted_otu_table.txt -o wf_jack/weighted_unifrac/sorted_otu_table_upgma.tre 
beta_diversity.py -i wf_jack/rarefaction/ -o wf_jack/weighted_unifrac/rare_dm/  -m weighted_unifrac  -t picked_otus_uclust/rep_set.tre
upgma_cluster.py -i wf_jack/weighted_unifrac/rare_dm/ -o wf_jack/weighted_unifrac/rare_upgma/ 
consensus_tree.py -i wf_jack/weighted_unifrac/rare_upgma/ -o wf_jack/weighted_unifrac/rare_upgma_consensus.tre 
tree_compare.py -s wf_jack/weighted_unifrac/rare_upgma/ -m wf_jack/weighted_unifrac/rare_upgma_consensus.tre -o wf_jack/weighted_unifrac/upgma_cmp/ 

mkdir wf_jack/unweighted_unifrac/
upgma_cluster.py -i wf_jack/unweighted_unifrac_sorted_otu_table.txt -o wf_jack/unweighted_unifrac/sorted_otu_table_upgma.tre 
beta_diversity.py -i wf_jack/rarefaction/ -o wf_jack/unweighted_unifrac/rare_dm/  -m unweighted_unifrac  -t picked_otus_uclust/rep_set.tre
upgma_cluster.py -i wf_jack/unweighted_unifrac/rare_dm/ -o wf_jack/unweighted_unifrac/rare_upgma/ 
consensus_tree.py -i wf_jack/unweighted_unifrac/rare_upgma/ -o wf_jack/unweighted_unifrac/rare_upgma_consensus.tre 
tree_compare.py -s wf_jack/unweighted_unifrac/rare_upgma/ -m wf_jack/unweighted_unifrac/rare_upgma_consensus.tre -o wf_jack/unweighted_unifrac/upgma_cmp/ 