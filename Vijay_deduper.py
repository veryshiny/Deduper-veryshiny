#!/usr/bin/env python

#command to run: 

#test file! WRITE THE STATEMENT
#./Vijay_deduper.py -f unit_test_folder/test_input_file.sam -u STL96.txt -o output.sam

import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import re

def get_args():
    parser = argparse.ArgumentParser(description="Reference Based PCR Duplicate Removal tool")
    parser.add_argument("-f", "--file", help="input SAM file", type= str)
    parser.add_argument("-u", "--umi", help="input UMI file", type= str)
    parser.add_argument("-o", "--outfile", help="output file", type= str)
    return parser.parse_args()

args = get_args()

UMI_file=open(args.umi,"r")
output= args.outfile
file_input=args.file

def retrieve_chr_num(line: str) -> int:
    '''Takes in a line of read feature data and outputs the chromosome number of the read'''
    chr_num= line.split('\t')[2]
    return chr_num

def retrieve_UMI(line: str) -> str:
    '''Takes in a line of read feature data and outputs the UMI of the read'''
    UMI = line.split('\t')[0].split(':')[-1]
    return UMI

def retrieve_strand(line: str) -> str:
    '''Takes in a line of read feature data and outputs the strand of the read'''
    flag=int(line.split('\t')[1])
            
    if ((flag & 16) != 16):
        strand="+"
    
    else:
        strand='-'
       
    return strand 

fake_cigar="NS500451:154:HWKTMBGXX:1:11101:94095:71756:AACGCCAT	0	2	76875967	36	1S72M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A"

def retrieve_position_clipping_cigar(line: str):
    '''Takes in a line of read feature data and outputs the actual left-most position of the read'''
    cigar=line.split('\t')[5]
    og_position=int(line.split('\t')[3])
    strand=retrieve_strand(line)
    dict_cigar={}
    dict_cigar['D']=0
    dict_cigar['N']=0 
    if strand=="+":
        if 'S' in cigar:
            if cigar.split('S')[0].isdigit():
                actual_position = og_position - int(cigar.split('S')[0]) 
                
            else:
                actual_position= og_position
        else:
            actual_position= og_position    
    else:
        match = re.findall(r'(\d+)(\w)', cigar)
        if match[-1][1]=='S':
            S_value = int(match [-1][0] )
        else:
            S_value = 0
        for tuple in match:
            if tuple[1] not in dict_cigar:
                dict_cigar[tuple[1]]=int(tuple[0])
            else:
                dict_cigar[tuple[1]]=dict_cigar[tuple[1]]+int(tuple[0])
        #print(dict_cigar)                   
        actual_position= og_position-1+S_value+dict_cigar["M"]+dict_cigar["D"]+dict_cigar["N"]
 
    return actual_position

#print(retrieve_position_clipping_cigar(fake_cigar))

set_UMIs=set()
for line in UMI_file:
    line= line.strip('\n')
    set_UMIs.add(line)

#print(set_UMIs)


set_umi_str_position=set()
chr_num= ""

with open(output,"w") as op: ##to get the multiple lines on a single line
    with open(file_input,"r") as file:
        for line in file:
            if line.startswith('@'):
                op.writelines(line)
            else:
                UMI=retrieve_UMI(line)
                strand=retrieve_strand(line)
                actual_position= retrieve_position_clipping_cigar(line)
                if UMI not in set_UMIs:
                    continue
                if chr_num != retrieve_chr_num(line):
                    chr_num = retrieve_chr_num(line)
                    set_umi_str_position=set()
                    set_umi_str_position.add((UMI,strand,actual_position))
                    op.writelines(line)
                else:
                    if (UMI,strand,actual_position) not in set_umi_str_position:
                        set_umi_str_position.add((UMI,strand,actual_position))
                        op.writelines(line)
                
