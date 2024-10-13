# **Pseudocode for Deduper** :rose: 

 A strategy for writing a Reference Based PCR Duplicate Removal tool.

 Given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). 
 
 A strategy that avoids loading everything into memory.

# **Problem**

Whet defines a PCR Duplicate:

- Same alignment position 
    • Chromosome            RNAME (SAM col 3)
    • 5' start of read          POS (SAM col 4) + CIGAR starting/ending with 'S' (SAM col 6)
    • Strand (strand specific?)     FLAG (SAM col 2 : BITFLAG 16)
- Same Unique Molecular Index (UMI or “randomer”)       QNAME (SAM col 1)

These PCR duplicates are artefacts of random extra amplification of certain reads during library prep. For downstream analyses they might cause problems while doing transcript abundance analysis while looking for genes with higher or lower expression.
 
# **Examples:** :circus_tent:
        Include a properly formated sorted input sam file
        Include a properly formated expected output sam file

## ***Case 1*** : **2 reads are not duplicates (are on different chromosomes)**

## ***Case 2*** : **2 reads are not duplicates (same chromosome but not same UMI)**

## ***Case 3*** : **2 reads are not duplicates (same chromosome, same UMI, not same strand)**

## ***Case 4*** : **2 reads are not duplicates (same chromosome, same UMI, same strand, not same position)**

## ***Case 5*** : **2 reads are not duplicates (same chromosome, same UMI, same strand, same position :: after soft-clipping they're not the same!)**

## ***Case 6*** : **2 reads are duplicates (same chromosome, same UMI, same strand, not same position:: but after soft-clipping they're the same!)**

## ***Case 7*** : **2 reads are duplicates**


# **Pseudocode** :round_pushpin:

## Step 1 : Sorting the files

We need to sort the files such that all the chromosomes are grouped together and then they're sorted by UMI.

Two ways we could do this:

- Use command line tools to sort by UMI and then chromosomes. 
- Use samtools to sort by chromosomes and then command-line sort by UMI (if first too slow)

## Step 2 : Python
 
![pseudocode](/home/varsheni/bgmp/bioinfo/Bi624/Deduper-veryshiny/pseudocode.png)

# **High-Level Functions**
    Determine high level functions
        Description
        Function headers
        Test examples (for individual functions)
        Return statement
