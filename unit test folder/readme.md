# Unit test summary statistics

|File| Number of reads |Number of Duplicate reads |
|---|---|
|[test_input_file.sam](test_input_file.sam)|14|4|
|[test_output_file.sam](test_output_file.sam)|12|2|

QNAME is the first column of the input and output files after the header lines containing @

### ***Case 1*** : **2 reads are not duplicates (are on different chromosomes)**

QNAME : 2_reads_are_on_different_chr

### ***Case 2*** : **2 reads are not duplicates (same chromosome but not same UMI)**

QNAME : 2_reads_are_on_same_chr_different_UMI

### ***Case 3*** : **2 reads are not duplicates (same chromosome, same UMI, not same strand)**

QNAME : 2_reads_are_on_same_chr_same_UMI_diff_strand

### ***Case 4*** : **2 reads are not duplicates (same chromosome, same UMI, same strand, not same position)**

QNAME: same_chr_same_UMI_same_strand_diff_position

### ***Case 5*** : **2 reads are not duplicates (same chromosome, same UMI, same strand, same position :: after soft-clipping they're not the same!)**

QNAME: same_chr_same_UMI_same_strand_diff_position_after_soft_clip

### ***Case 6*** : **2 reads are duplicates (same chromosome, same UMI, same strand, not same position:: but after soft-clipping they're the same!)**

QNAME: dup_same_chr_same_UMI_same_strand_same_position_after_soft_clip

### ***Case 7*** : **2 reads are duplicates**

QNAME: dup_same_chr_same_UMI_same_strand_same_position

### ***Case 8*** : **read has unknown UMI**

QNAME: unknown_UMI