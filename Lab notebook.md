# Lab Notebook :cherry_blossom:

# Day 1

## initial problem analysis

maybe running through samtools sort and then sorting by umi per chromosome would work 

```bash
samtools sort  -O sam -o output.sam  test.sam
grep -v "^@" output.sam | sort -k3,1  -t ':'  -k 8,2 > output1.sam
```

Its wrong it sorts by UMI first somehow.

```bash
grep -v "^@" test.sam | sort -t ':'  -k 8 | sort -k 3 -s > output2.sam
```

# CIGAR

|Letter|Meaning|Consume query |Consume reference|5' |3'|
|--|--|--|--|--|--|
|M | match and mismatch|Yes|Yes|||
|S |soft clip|Yes|No|The beginning|The ending|
|I |insertion|Yes|No|||
|D| deletion|No|Yes|||
|N|Alignment gap|No|Yes|||


# Pairwise, choice of dupe and error correction.