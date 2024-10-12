# **Pseudocode for Deduper** :rose: :sunflower:

 A strategy for writing a Reference Based PCR Duplicate Removal tool.

 Given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory.

# **Problem**

Whet defines a PCR Duplicate:

Same alignment position 
    • Chromosome             RNAME (SAM col 3)
    • 5' start of read       POS (SAM col 4) + CIGAR (SAM col 6)
    • Strand (strand specific?) FLAG (SAM col 2 : BITFLAG 16)
Same Unique Molecular Index (UMI or “randomer”)  QNAME (SAM col 1)

 
# **Examples:** :circus_tent:
        Include a properly formated sorted input sam file
        Include a properly formated expected output sam file
    
# **Pseudocode**
    Develop your algorithm using pseudocode

# **High-Level Functions**
    Determine high level functions
        Description
        Function headers
        Test examples (for individual functions)
        Return statement
