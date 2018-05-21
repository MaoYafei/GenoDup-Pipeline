# Genome-duplication-calculation
  This python script can be used to detect Whole-genome duplication (WGD) in the dS based method.

# Run the program
  1. Caculate dS values based on gene family data (Paralogs dS values)

  python main.py -s Nuclear_sequence_file -p Protein_sequence_file -g Gene_family_file -n Maximum_number_of_gene_family

  2. Caculate dS values based on gene pair data (Anchor dS values)

  python main.py -s Nuclear_sequence_file -p Protein_sequence_file -c Gene_pair_file

# Input
  1. Nuclear_sequence_file: it contains all the nuclear sequences in your analysis (fasta format).
    eg:
        >gene1
          ATCG
        >gene2
          ATCC
        ...
  2. Protein_sequence_file: it contains all the protein sequences in your analysis (fasta format).
    eg:
        >gene1
          PAPA
        >gene2
          PAPA
        ...
  3. Gene_family_file: it contains the gene family cluster (usually be produced by OrthoMCL).
     eg:
        led1: gene1,gene2,gene3
        led2: gene3,gene4
        ...
  4. Gene_pair_file: it contains two Ohnologs in two colums separated by tab (could be produced by MCScanX or OrthoMCL or i-ADHoRe).
     eg:
        gene1 gene2
        gene3 gene4
        ...
  5. Maximum_number_of_gene_family: Maximum number of gene family which you want to analyze, only use with -g (suggest: 5-10)

# Output

  1. pairwise directory: including all gene pairs sequence (fasta format).
  2. PAML_result: including all codeml output files.

# Visualization

  Running Geno_Dup_Cal_plot.sh
  
# Citation


# Others
  Please read below literatures for basics on dS value caculation for WGD detectation:
    1.
    2.
    3.
    4.

