# **Pseudocode for Deduper** :rose: 

 A strategy for writing a Reference Based PCR Duplicate Removal tool.

 Given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). 
 
 A strategy that avoids loading everything into memory.

# **Problem**

Whet defines a PCR Duplicate:

- **Same alignment position** 

| Identifier | Description |
| --- | --- |
| Chromosome   | RNAME (SAM col 3)  |
| 5' start of read  | POS (SAM col 4) + CIGAR  starting/ending with 'S' (SAM col 6)  |
| Strand | FLAG (SAM col 2 : BITFLAG 16) |
       
- **Same Unique Molecular Index**       
QNAME (SAM col 1)



PCR duplicates are identified by being present on the same position on the same chromosome AND having the same UMI, which is an index randomly inserted during library prep. The UMI is a way to tell whether the read was actually a technical replicate or a biological one, with more reassurance. If the reads mapped to same positions but they do not have the same UMIs, they could be biological replicates and indicate actual increase in gene expression. 

Duplicates are good during alignment as it increases coverage of particular regions. This is why its more important to remove it after aligning, also removing them during alignment could be very time and memory consuming as there would need to be 1-1 comparisions of millions of reads. 

Clipping can help improve the quality of alignments by focusing on the high-confidence, well-aligned portions of the reads and discarding or ignoring the low-quality or non-aligning parts. Soft clipping retains the unaligned portions of the read in the alignment file but does not use them in the actual alignment. This means that the unaligned bases are not considered when calculating alignment scores or making downstream analyses. This is accounted for in the CIGAR string indicated by an S. We need to make sure that the actual positions of 2 reads with the same UMI, same chromosome and same strand have the same position by accounting for the CIGAR value as well.

These PCR duplicates are artefacts of random extra amplification of certain reads during library prep. For downstream analyses they might cause problems while doing transcript abundance analysis while looking for genes with higher or lower expression.
 
# **Examples:** :circus_tent:
        
Formatted sorted test input sam file
Formatted sorted expected output sam file

### ***Case 1*** : **2 reads are not duplicates (are on different chromosomes)**

### ***Case 2*** : **2 reads are not duplicates (same chromosome but not same UMI)**

### ***Case 3*** : **2 reads are not duplicates (same chromosome, same UMI, not same strand)**

### ***Case 4*** : **2 reads are not duplicates (same chromosome, same UMI, same strand, not same position)**

### ***Case 5*** : **2 reads are not duplicates (same chromosome, same UMI, same strand, same position :: after soft-clipping they're not the same!)**

### ***Case 6*** : **2 reads are duplicates (same chromosome, same UMI, same strand, not same position:: but after soft-clipping they're the same!)**

### ***Case 7*** : **2 reads are duplicates**

# **Pseudocode** :round_pushpin:

## Step 1 : Sorting the files

We need to sort the files such that all the chromosomes are grouped together and then they're sorted by UMI.

Two ways we could do this:

- Use command line tools to sort by UMI and then chromosomes. 
- Use samtools to sort by chromosomes and then command-line sort by UMI (if first too slow)

## Step 2 : Python
 
![pseudocode](pseudocode.png)

# **High-Level Functions**
    Determine high level functions
        Description
        Function headers
        Test examples (for individual functions)
        Return statement

## *Function 1: Chromosome number calculator*
