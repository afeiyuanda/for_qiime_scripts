#!/usr/bin/python 
#-*-coding:utf-8-*-
'''
Author:	Liu wei
Data:	2014.2.12
Usage:	python blast_fomat8_file_dereplication total.out statistic.txt
delet the replication lines by the first column,keep the biggest identity line, then calculate the second column's names
total.out files:
HWI-ST1330:42:C1CW6ACXX:5:1311:17247:59030      AFE_2629        86.81   91      12      0       1       91      91      1       9e-20   85.7
HWI-ST1330:42:C1CW6ACXX:5:2202:6598:76375       AFE_2627        85.11   94      14      0       4       97      614     707     9e-17   75.8
HWI-ST1330:42:C1CW6ACXX:5:2203:12234:9454       AFE_3142        89.80   49      5       0       3       51      350     302     2e-11   58.0
HWI-ST1330:42:C1CW6ACXX:5:1111:15298:74233      AFE_2618        82.00   100     18      0       1       100     440     539     8e-11   56.0
HWI-ST1330:42:C1CW6ACXX:5:2114:5283:44917       AFE_2624        90.00   90      9       0       1       90      307     396     2e-26    107
HWI-ST1330:42:C1CW6ACXX:5:2316:13351:81878      AFE_2550        92.86   42      3       0       4       45      847     888     5e-12   60.0
HWI-ST1330:42:C1CW6ACXX:5:1202:10227:85059      AFE_2550        86.90   84      11      0       1       84      69      152     6e-18   79.8
HWI-ST1330:42:C1CW6ACXX:5:2110:20835:42524      AFE_3144        85.54   83      12      0       4       86      196     114     5e-15   69.9
HWI-ST1330:42:C1CW6ACXX:5:1213:4311:46673       AFE_3209        88.41   69      8       0       1       69      231     299     3e-16   73.8
HWI-ST1330:42:C1CW6ACXX:5:2105:16808:71544      AFE_2617        83.87   93      15      0       1       93      1107    1015    8e-14   65.9
HWI-ST1330:42:C1CW6ACXX:5:2108:18540:65481      AFE_2554        82.29   96      17      0       4       99      15      110     8e-11   56.0
'''

import sys,os,subprocess

def usage():
    print """
			Usage:	python blast_fomat8_file_dereplication total.out statistic.txt
			delet the replication lines by the first column, then calculate the second column's names
			total.out files:
			HWI-ST1330:42:C1CW6ACXX:5:1311:17247:59030      AFE_2629        86.81   91      12      0       1       91      91      1       9e-20   85.7
			HWI-ST1330:42:C1CW6ACXX:5:2202:6598:76375       AFE_2627        85.11   94      14      0       4       97      614     707     9e-17   75.8
			HWI-ST1330:42:C1CW6ACXX:5:2203:12234:9454       AFE_3142        89.80   49      5       0       3       51      350     302     2e-11   58.0
			
			total.out.titles_no_replication files:
			1 HWI-ST1330:42:C1CW6ACXX:5:1101:10081:22438
			2 HWI-ST1330:42:C1CW6ACXX:5:1101:10084:56075
			1 HWI-ST1330:42:C1CW6ACXX:5:1101:10087:53592
			2 HWI-ST1330:42:C1CW6ACXX:5:1101:10101:20440
		"""

def identity_biggest_line( lines ):
	max = 0
	for i in lines:
		#use bit score to identity goodest blast hit
		identity = i.split( '\t' )[-3]
		if float(identity) > float(max):
			max = float(identity)
			target_line = i
	return target_line
	
def __main__():

	try:
		blast_output = sys.argv[1]
	except:
		usage()
		sys.exit(1)
	
	f = open( blast_output )
	blast_lines = f.readlines()
	f.close()
	query_no_replication = blast_output + "_deduplicate"
	output = open( query_no_replication, 'w' )
	temp_lines_list = []
	query_list = []
	title_line = blast_lines[0]
	output.write(title_line)
	for line in blast_lines[1:]:
		if line.split('\t')[0] not in query_list:
			if temp_lines_list != []:
				output.write(identity_biggest_line(temp_lines_list))
				temp_lines_list = []
			query_list.append(line.split('\t')[0])
			temp_lines_list.append(line)
		else:
			temp_lines_list.append(line)

	output.close()
if __name__ == "__main__": __main__()