# Deduper Super Deluxe!

Given a SAM file of uniquely mapped reads, and a text file containing the known UMIs OR using random-mers as the UMIs instead; remove all PCR duplicates (retain only a single copy of each read). 

The Python code assumes a sorted sam file (you can use `samtools sort` first!)
- The script called `Vijay_deduper_super_deluxe.py` accounts for: 
    - all possible CIGAR strings (including adjusting for soft clipping, etc.)
    - Strand
    - Single-end reads
    - Known UMIs vs randomers 
    - Error correction of known UMIs 
    - Highest quality and longest duplicate written to file 
- Considerations:
    - Millions of reads – avoids loading everything into memory! 
    - Runs in less than 2 minutes (even with all the new add-ons)
    - Outputs a properly formatted SAM file

- **NEW FEATURES COMING SOON**
    - Single-end and paired-end 

    
- Includes the following argparse options
    - ```-f```, ```--file```: designates absolute file path to sorted sam file
    - ```-o```, ```--outfile```: designates absolute file path to deduplicated sam file
    - ```-u```, ```--umi```: designates file containing the list of UMIs
    - ```-e```, ```--umierrorcorrection```: designates whether the user wants error correction of UMIs
    - ```-r```, ```--randomers```: designates that the UMIs should be considered as randomers of a certain length instead of a provided list
    - ```-p```, ```--paired_end```: designates absolute file path to sorted sam file with paired-end data
    - ```-h```, ```--help```: prints a USEFUL help message
        - That is, the code is able to run (in a single step) if given a command in the format:
          ```
          ./Vijay_deduper_super_deluxe.py -u STL96.txt -f <in.sam> -o <out.sam>
          ```



