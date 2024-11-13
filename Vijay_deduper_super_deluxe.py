#!/usr/bin/env python

#command to run: 

#./Vijay_deduper_super_deluxe.py -f unit_test_folder/test_input_file.sam -u STL96.txt -o output.sam
#Remember that the sam file input needs to be sorted by chromosome!
#python -m cProfile -o output_filename.pstats ./Vijay_deduper_super_deluxe.py -u STL96.txt -f sorted_C1_SE_uniqueAlign.sam -o deduper_C1_SE_uniqAlign.sam

import argparse
import re
import itertools

def get_args():
    parser = argparse.ArgumentParser(description="Reference Based PCR Duplicate Removal tool")
    parser.add_argument("-f", "--file", help="input SAM file", type= str)
    parser.add_argument("-u", "--umi", help="input UMI file", type= str)
    parser.add_argument("-o", "--outfile", help="output file", type= str)
    parser.add_argument("-p", "--paired_end", help="sam file for paired end data", action='store_true')
    parser.add_argument("-e", "--umierrorcorrection", help="use flag if you want umis to be error corrected; uses average phred score and length of reads as filter", action='store_true')
    parser.add_argument("-r", "--randomers", help="use flag if you want to use randomers instead of a known UMI list", action='store_true')

    return parser.parse_args()

args = get_args()
UMI_error_correction= args.umierrorcorrection

paired_end=args.paired_end
UMI_file=open(args.umi,"r")
output= args.outfile
file_input=args.file
randomners=args.randomers

#def retrieve_chr_num(line: str) -> str:
#    '''Takes in a line of read feature data and outputs the chromosome number of the read'''
#    chr_num= line.split('\t')[2]
#    return chr_num

#def retrieve_UMI(line: str) -> str:
#    '''Takes in a line of read feature data and outputs the UMI of the read'''
#    UMI = line.split('\t')[0].split(':')[-1]
#    return UMI

def retrieve_strand(flag: int) -> str:
    '''Takes in a line of read feature data and outputs the strand of the read'''
    #flag=int(line.split('\t')[1])
            
    if ((flag & 16) != 16):
        strand="+"
    else:
        strand='-'
    return strand 

fake_cigar="NS500451:154:HWKTMBGXX:1:11101:94095:71756:AACGCCAT	0	2	76875967	36	1S72M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A"

pattern= re.compile(r'(\d+)(\w)')

def retrieve_position_clipping_cigar(cigar,og_position,strand: str):
    '''Takes in a line of read feature data and outputs the actual left-most position of the read'''
    #cigar=line.split('\t')[5]
    #og_position=int(line.split('\t')[3])
    S_value=0
    if strand=="+":
        if 'S' in cigar:
            number_S=cigar.split('S')[0]
            if number_S.isdigit():
                actual_position = og_position - int(number_S)                 
            else:
                actual_position= og_position
        else:
            actual_position= og_position    
    else:
        dict_cigar={'D':0,'N':0,'I':0,'S':0,'M':0,'X':0,
                '=':0,'H':0,'P':0}    
        match = pattern.findall(cigar)
        if cigar[-1]=='S':
            S_value = int(match [-1][0] )
        for tuple in match:
            dict_cigar[tuple[1]]=dict_cigar[tuple[1]]+int(tuple[0])
        #print(dict_cigar)                   
        actual_position= og_position-1+S_value+dict_cigar["M"]+dict_cigar["D"]+dict_cigar["N"] 
    return actual_position

#print(retrieve_position_clipping_cigar(fake_cigar,retrieve_strand(fake_cigar)))

def quality_finder(phred):
    sum_quality=0
    for letter in phred:
        sum_quality+=ord(letter)
    average_quality=sum_quality/len(phred)
    return average_quality


set_UMIs=set()
list_UMIs=[]
if randomners==False:
    for line in UMI_file:
        line= line.strip('\n')
        set_UMIs.add(line)
        list_UMIs.append(line)
else:
    Base_array="ATCG"
    permut= list(itertools.product(['A', 'C', 'T','G'], repeat=8))
    def convertTuple(tup):
            # initialize an empty string
        str = ''
        for item in tup:
            str = str + item
        return str
    
    for tuple in permut:
        set_UMIs.add(convertTuple(tuple))

#print(len(set_UMIs))

#print(set_UMIs)

#n_cols, n_rows = len(list_UMIs[0]), len(set_UMIs)
#matrix_for_correction = np.zeros((n_rows, n_cols), dtype=np.int32)
#df = pd.DataFrame(matrix_for_correction, index=list_UMIs)
#print(matrix_for_correction)

def error_correction(UMI: str):
    break_up_UMI=list(UMI)
   
    list_scores=[]
    for good_UMI in list_UMIs:
        score = 0
        for number in range(0,len(UMI)):
         if break_up_UMI[number]==good_UMI[number]:
             score +=1
        list_scores.append(score)
    if list_scores.count(max(list_scores))==1  and max(list_scores)>len(UMI)/2+2: 
        UMI_value_return= list_UMIs[list_scores.index(max(list_scores)) ]
    else:
        UMI_value_return=False
    return UMI_value_return

#print(error_correction("AATCTGTG"))

set_umi_position_positive=set()
set_umi_position_negative=set()

chr_num= ""
print_number_chromosome=0
dictionary_for_length_and_quality_positive={}
dictionary_for_length_and_quality_negative={}


if paired_end==False :
    with open(output,"w") as op: ##to get the multiple lines on a single line
        with open(file_input,"r") as file:
            number_unknown_UMI= 0
            number_dupes=0
            number_unique=0
            number_header=0
            count_line_number=0
            corrected_UMIs=0
            #header=True
            for line in file:
                
                if line.startswith('@'):
                    op.write(line)
                    number_header+=1
                else:
                    #header = False
                    split_line=line.split('\t')
                    UMI=split_line[0].split(':')[-1]
                    if chr_num != split_line[2]:
                        
                        
                        for key in dictionary_for_length_and_quality_positive.keys():
                            line_to_write=dictionary_for_length_and_quality_positive[key]
                            #print(line_to_write)
                            op.write(line_to_write[2])
                        for key in dictionary_for_length_and_quality_negative.keys():
                            line_to_write=dictionary_for_length_and_quality_negative[key]
                            op.write(line_to_write[2])    
                        set_umi_position_negative=set()
                        set_umi_position_positive=set()
                        dictionary_for_length_and_quality_positive={}
                        dictionary_for_length_and_quality_negative={}
                        print(chr_num,print_number_chromosome,sep="\t")
                        chr_num=split_line[2]
                        print_number_chromosome=0
                    if UMI in set_UMIs:
                        strand=retrieve_strand(int(split_line[1]))
                        
                        actual_position= retrieve_position_clipping_cigar(split_line[5],int(split_line[3]),strand)
                        quality=quality_finder(split_line[10])
                        if strand=='+': 
                            if(UMI,actual_position) not in set_umi_position_positive:
                                
                                set_umi_position_positive.add((UMI,actual_position))
                                dictionary_for_length_and_quality_positive[(UMI,actual_position)] = (len(split_line[9]),quality,line) # length, quality,line
                                number_unique+=1
                                print_number_chromosome+=1
                            else:
                                number_dupes+=1
                                if dictionary_for_length_and_quality_positive[(UMI,actual_position)][1] < quality :
                                    dictionary_for_length_and_quality_positive[(UMI,actual_position)]= (len(split_line[9]),quality,line)
                                elif dictionary_for_length_and_quality_positive[(UMI,actual_position)][1] == quality and dictionary_for_length_and_quality_positive[(UMI,actual_position)][0]<len(split_line[9]):
                                    dictionary_for_length_and_quality_positive[(UMI,actual_position)]= (len(split_line[9]),quality,line)


                        elif (UMI,actual_position) not in set_umi_position_negative:
                            #op.write(line)
                           
                            set_umi_position_negative.add((UMI,actual_position))
                            dictionary_for_length_and_quality_negative[(UMI,actual_position)] = (len(split_line[9]),quality,line) # length, quality,line
                            
                            number_unique+=1
                            print_number_chromosome+=1
                        else:
                            number_dupes+=1
                            
                            if dictionary_for_length_and_quality_negative[(UMI,actual_position)][1] < quality :
                                dictionary_for_length_and_quality_negative[(UMI,actual_position)]= (len(split_line[9]),quality,line)
                            elif dictionary_for_length_and_quality_negative[(UMI,actual_position)][1] == quality and dictionary_for_length_and_quality_negative[(UMI,actual_position)][0]<len(split_line[9]):
                                dictionary_for_length_and_quality_negative[(UMI,actual_position)]= (len(split_line[9]),quality,line)
                    elif UMI_error_correction==True:
                        UMI = error_correction(UMI)
                        #print(UMI)
                        if UMI== False:
                            number_unknown_UMI+=1
                        else:
                            corrected_UMIs+=1
                            strand=retrieve_strand(int(split_line[1]))
                            actual_position= retrieve_position_clipping_cigar(split_line[5],int(split_line[3]),strand)
                            quality=quality_finder(split_line[10])
                            if strand=='+': 
                                if(UMI,actual_position) not in set_umi_position_positive:
                                    #op.write(line)
                                    set_umi_position_positive.add((UMI,actual_position))
                                    dictionary_for_length_and_quality_positive[(UMI,actual_position)] = (len(split_line[9]),quality,line) # length, quality,line
                                    number_unique+=1
                                    print_number_chromosome+=1
                                else:
                                    number_dupes+=1
                                    if dictionary_for_length_and_quality_positive[(UMI,actual_position)][1] < quality :
                                        dictionary_for_length_and_quality_positive[(UMI,actual_position)]= (len(split_line[9]),quality,line)
                                    elif dictionary_for_length_and_quality_positive[(UMI,actual_position)][1] == quality and dictionary_for_length_and_quality_positive[(UMI,actual_position)][0]<len(split_line[9]):
                                        dictionary_for_length_and_quality_positive[(UMI,actual_position)]= (len(split_line[9]),quality,line)


                            elif (UMI,actual_position) not in set_umi_position_negative:
                                #op.write(line)
                                set_umi_position_negative.add((UMI,actual_position))
                                dictionary_for_length_and_quality_negative[(UMI,actual_position)] = (len(split_line[9]),quality,line) # length, quality,line
                                number_unique+=1
                                print_number_chromosome+=1
                            else:
                                number_dupes+=1
                                if dictionary_for_length_and_quality_negative[(UMI,actual_position)][1] < quality :
                                    dictionary_for_length_and_quality_negative[(UMI,actual_position)]= (len(split_line[9]),quality,line)
                                elif dictionary_for_length_and_quality_negative[(UMI,actual_position)][1] == quality and dictionary_for_length_and_quality_negative[(UMI,actual_position)][0]<len(split_line[9]):
                                    dictionary_for_length_and_quality_negative[(UMI,actual_position)]= (len(split_line[9]),quality,line)
                                    


                    else:
                         number_dupes+=1

        for key in dictionary_for_length_and_quality_positive.keys():
            line_to_write=dictionary_for_length_and_quality_positive[key]
            op.write(line_to_write[2])
        for key in dictionary_for_length_and_quality_negative.keys():
            line_to_write=dictionary_for_length_and_quality_negative[key]
            op.write(line_to_write[2]) 

print(chr_num,print_number_chromosome,sep="\t")
print("\nHeader lines =",number_header,"\nUnique reads =",number_unique,"\nNumber of unknown UMIs =",number_unknown_UMI,"\nNumber of duplicates =",number_dupes,"\nCorrected UMIs =",corrected_UMIs)      
            