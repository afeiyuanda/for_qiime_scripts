#-*-coding:utf-8-*-
'''
Author:    Liu wei
Data:    2018.01.15
'''
#need cogent/ and fastq.py
import os,sys,re,argparse
from fastq import MinimalFastqParser

#Reverse complementary sequence
def reverse_com_seq( seq ):
    new_seq = ''.join(["ATCGNFXRYMKSWHBVD"["TAGCNFXRYMKSWHBVD".index(n)] for n in seq[::-1]])
    return new_seq

def __main__():
    try:
        input_seqs_ul = sys.argv[1]
        target_seq = sys.argv[2]
        output_dir = sys.argv[3]
    except:
        sys.exit(1)
    
    outfile_name = os.path.basename(input_seqs_ul).split('.')[0]
    output_fq =  output_dir + '/'+outfile_name
    outf = open( output_fq , 'w')
    
    greped_num = 0
    for seq_id, seq, qual in MinimalFastqParser(input_seqs_ul, strict=False):
        if target_seq in seq:
            outf.write('@%s\n%s\n+\n%s\n' % (seq_id, seq, qual))
            greped_num += 1
        elif reverse_com_seq(target_seq) in seq:
            outf.write('@%s\n%s\n+\n%s\n' % (seq_id+'|reverse', reverse_com_seq(seq), qual))
            greped_num += 1
        else:
            pass
    print 'greped_num is %s' % greped_num
    outf.close()

if __name__ == "__main__" : __main__()