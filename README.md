## ProteinMPNN_utils
**Various utility functions for reading ProteinMPNN-formatted FASTA files and outputting regular FASTA files (e.g. for AF2), storing all the information in pandas dataframes, and doing a few operations.**
###### Robert Alberstein, 2022

**FUNCTIONS:**
- MPNN_fasta_to_df(fname):
  - reads in the MPNN-formatted FASTA file fname and returns a pandas DF
- fasta_to_df(fname):
  - reads in a standard format FASTA file fname and returns a pandas DF
- write_df_to_fasta(df, fname):
  - writes the DataFrame df to FASTA file fname
- add_loop(seq, loop_seq):
  - replaces the "/" in a sequence seq with the string loop_seq
- calc_num_residue_df(df, res_1letter):
  - given a dataframe and a 1-letter residue, count the number of instances of that residue in all sequences and add this as a new column to the DF
- calc_seq_identity_array(seq1, seq2):
  - calculate a binary array of position identities between two sequences
