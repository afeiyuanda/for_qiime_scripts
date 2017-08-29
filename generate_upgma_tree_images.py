#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2015.1.4
Usage:	python generate_upgma_tree_images.py
'''

import sys,subprocess,os

def usage():
    print """
			Usage:	python generate_upgma_tree_images.py
					Firstly, we change the tree format to NEXUS format;
					Secondly, we use TreeGraph_2.3.0 generate eps images;
					Finally, we convert eps to png.
		"""

'''
#NEXUS
begin trees;
tree tree1 = [&r] ((DS.NT1:0.130508269921, DS.NT2:0.130508269921):0.0954473085839, ((DS.SD1:0.0534318208235, DS.YN2:0.0534318208235):0.0407428494751, (DS.SD2:0.064180676306, DS.YN1:0.064180676306):0.0299939939926):0.131780908207);
end;
'''
		
def format2NEXUS_then_plot( ul ):
	( file_path, file_name ) = os.path.split( ul )
	file_name_prefix = os.path.splitext( file_name )[0]
	new_ul = file_path + '/' + file_name_prefix + '_nexus.tre'
	tre_line = open( ul ).readlines()[0]
	output = open( new_ul, 'w' )
	output.writelines( '#NEXUS\nbegin trees;\n' )
	output.writelines( 'tree tree1 = [&r] ' + tre_line + 'end;' )
	output.close()
	print "Tree format changed successfully!"
	
	subprocess.call( "tgf -t " + new_ul + ' \proof false', stdout=subprocess.PIPE, shell=True )
	subprocess.call( "tgf -p " + file_path + '/' + file_name_prefix + '_nexus.tgf', stdout=subprocess.PIPE, shell=True )
	#Linux command: convert ds_upgma_NEXUS.eps ds_upgma_NEXUS.png
	subprocess.call( "convert " + file_path + '/' + file_name_prefix + '_nexus.eps ' + file_path + '/' + file_name_prefix + '_nexus.png', stdout=subprocess.PIPE, shell=True )
	print "Tree image generated successfully!"
	
def __main__():
	usage()
	indir = sys.argv[1]
	ul1 = indir+'/wf_jack/weighted_unifrac/upgma_cmp/'
	ul2 = indir+'/wf_jack/unweighted_unifrac/upgma_cmp/'
	new_ul1 = indir+'/UPGMA/weighted_unifrac_upgma/'
	new_ul2 = indir+'/UPGMA/unweighted_unifrac_upgma/'
	os.system( 'mkdir -p ' + new_ul1 )
	os.system( 'mkdir -p ' + new_ul2 )
	
	os.system( 'cp ' + ul1 + 'master_tree.tre ' + new_ul1 )
	os.system( 'cp ' + ul1 + 'jackknife_named_nodes.tre ' + new_ul1 )
	
	os.system( 'cp ' + ul2 + 'master_tree.tre ' + new_ul2 )
	os.system( 'cp ' + ul2 + 'jackknife_named_nodes.tre ' + new_ul2 )
	
	format2NEXUS_then_plot( ul1 + '/master_tree.tre' )
	format2NEXUS_then_plot( ul1 + '/jackknife_named_nodes.tre' )
	format2NEXUS_then_plot( ul2 + '/master_tree.tre' )
	format2NEXUS_then_plot( ul2 + '/jackknife_named_nodes.tre' )
	
	os.system( 'cp ' + ul1 + '*.png ' + new_ul1 )
	os.system( 'cp ' + ul2 + '*.png ' + new_ul2 )
	
if __name__ == "__main__": __main__()