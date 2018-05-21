# Genome-duplication-calculation
  This python script can be used to detect Whole-genome duplication (WGD) in the dS based method.

# Run the program
  1. Caculate dS values based on gene family data (Paralogs dS values)

  python main.py -s Nuclear_sequence_file -p Protein_sequence_file -g Gene_family_file -n Maximum_number_of_gene_family

  2. Caculate dS values based on gene pair data (Anchor dS values)

  python main.py -s Nuclear_sequence_file -p Protein_sequence_file -c Gene_pair_file

# Output

  1. pairwise directory: including all gene pairs sequence (fasta format).
  2. PAML_result: including all codeml output files.

# Visualization

  Running Geno_Dup_Cal_plot.sh

