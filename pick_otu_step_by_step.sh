#!/bin/bash
# interactive commands are commented out

# Pick Reference OTUs (prefilter) command 
parallel_pick_otus_uclust_ref.py -i split_library_output/seqs.fna -o picked_otus_uclust/prefilter_otus/ -r /home/liuwei/qiime_software/gg_otus-12_10-release/rep_set/97_otus.fasta -T --jobs_to_start 4 --similarity 0.6 --enable_rev_strand_match

# Filter prefilter failures from input command 
filter_fasta.py -f split_library_output/seqs.fna -o picked_otus_uclust/prefilter_otus//prefiltered_seqs.fna -s picked_otus_uclust/prefilter_otus//seqs_failures.txt ¨Cn

# Pick Reference OTUs command 
parallel_pick_otus_uclust_ref.py -i picked_otus_uclust/prefilter_otus//prefiltered_seqs.fna -o picked_otus_uclust/step1_otus -r /home/liuwei/qiime_software/gg_otus-12_10-release/rep_set/97_otus.fasta -T --jobs_to_start 4 --similarity 0.97 --enable_rev_strand_match

# Generate full failures fasta file command 
filter_fasta.py -f picked_otus_uclust/prefilter_otus//prefiltered_seqs.fna -s picked_otus_uclust/step1_otus/prefiltered_seqs_failures.txt -o picked_otus_uclust/step1_otus/failures.fasta

# Pick rep set command 
pick_rep_set.py -i picked_otus_uclust/step1_otus/prefiltered_seqs_otus.txt -o picked_otus_uclust/step1_otus/step1_rep_set.fna -f picked_otus_uclust/prefilter_otus//prefiltered_seqs.fna

# Pick de novo OTUs for new clusters command 
pick_otus.py -i picked_otus_uclust/step1_otus/subsampled_failures.fasta -o picked_otus_uclust/step2_otus/ -m uclust  --uclust_otu_id_prefix New.ReferenceOTU --similarity 0.97 --enable_rev_strand_match

# Pick representative set for subsampled failures command 
pick_rep_set.py -i picked_otus_uclust/step2_otus//subsampled_failures_otus.txt -o picked_otus_uclust/step2_otus//step2_rep_set.fna -f picked_otus_uclust/step1_otus/subsampled_failures.fasta

# Pick reference OTUs using de novo rep set command 
parallel_pick_otus_uclust_ref.py -i picked_otus_uclust/step1_otus/failures.fasta -o picked_otus_uclust/step3_otus/ -r picked_otus_uclust/step2_otus//step2_rep_set.fna -T --jobs_to_start 4 --similarity 0.97 --enable_rev_strand_match

# Create fasta file of step3 failures command 
filter_fasta.py -f picked_otus_uclust/step1_otus/failures.fasta -s picked_otus_uclust/step3_otus//failures_failures.txt -o picked_otus_uclust/step3_otus//failures_failures.fasta

# Pick de novo OTUs on step3 failures command 
pick_otus.py -i picked_otus_uclust/step3_otus//failures_failures.fasta -o picked_otus_uclust/step4_otus/ -m uclust  --uclust_otu_id_prefix New.CleanUp.ReferenceOTU --similarity 0.97 --enable_rev_strand_match

# Merge OTU maps command 
cat picked_otus_uclust/step1_otus/prefiltered_seqs_otus.txt picked_otus_uclust/step3_otus//failures_otus.txt picked_otus_uclust/step4_otus//failures_failures_otus.txt >> picked_otus_uclust/final_otu_map.txt

# Pick representative set for subsampled failures command 
pick_rep_set.py -i picked_otus_uclust/step4_otus//failures_failures_otus.txt -o picked_otus_uclust/step4_otus//step4_rep_set.fna -f picked_otus_uclust/step3_otus//failures_failures.fasta

# Make the otu table command 
make_otu_table.py -i picked_otus_uclust/final_otu_map_mc2.txt -o picked_otus_uclust/otu_table_mc2.biom

# Assign taxonomy command 
#parallel_assign_taxonomy_rdp.py -i picked_otus_uclust/rep_set.fna -o picked_otus_uclust/rdp_assigned_taxonomy -T --jobs_to_start 4 
parallel_assign_taxonomy_blast.py -i g_Phascolarctobacterium.fna  -e 0.00001 -O 8 -o assign_taxa_blast
assign_taxonomy.py -i g__Staphylococcus.fna -r $reference_seqs -t $reference_tax  -m blast
assign_taxonomy.py -o rdp_silva/ -i rep_set.fna  -r /home/liuwei/database/Silva_111_post/rep_set/97_Silva_111_rep_set.fasta  -t /home/liuwei/database/Silva_111_post/taxonomy/97_Silva_111_taxa_map_RDP_6_levels.txt -m rdp

# Add taxa to OTU table command by qiime 1.7
add_metadata.py -i picked_otus_uclust/otu_table_mc2.biom --observation_mapping_fp picked_otus_uclust/rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o picked_otus_uclust/otu_table_mc2_w_tax.biom --sc_separated taxonomy --observation_header OTUID,taxonomy

# Add taxa to OTU table command by qiime 1.8
biom add-metadata -i picked_otus_uclust/otu_table_mc2.biom --observation-metadata-fp  picked_otus_uclust/rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o picked_otus_uclust/otu_table_mc2_w_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy

# Align sequences command 
parallel_align_seqs_pynast.py -i picked_otus_uclust/rep_set.fna -o picked_otus_uclust/pynast_aligned_seqs -T --jobs_to_start 4 

# Filter alignment command 
filter_alignment.py -o picked_otus_uclust/pynast_aligned_seqs -i picked_otus_uclust/pynast_aligned_seqs/rep_set_aligned.fasta 

#Filter reads in the biom file which can not alligned 
filter_otus_from_otu_table.py -i picked_otus_uclust/otu_table_mc2_w_tax.biom -e picked_otus_uclust/pynast_aligned_seqs/rep_set_failures.fasta -o otu_table_mc2_w_tax_no_pynast_failures.biom

# Build phylogenetic tree command 
make_phylogeny.py -i picked_otus_uclust/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta -o picked_otus_uclust/rep_set.tre 
