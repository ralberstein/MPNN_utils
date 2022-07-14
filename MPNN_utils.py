"""
Various utility functions for reading MPNN-formatted FASTA files and outputting regular FASTA files (e.g. for AF2), storing all the information in pandas dataframes, and doing a few operations.
Robert Alberstein, 2022
"""
import pandas as pd

##### FUNCTIONS #####
def MPNN_fasta_to_df(fname):
	"""
	Reads in an MPNN output fasta file and returns a pandas df with columns [T, sample, score, seq_recovery, sequence]
	"""
	all_data = []
	with open(fname, "r") as fin:
		# Store info for original sequence
		orig_info = fin.readline()[1:-1]  # omit starting '>' and newline character
		orig_score = orig_info[1:-1].split()[1].split("=")[1][:-1]
		orig_seq = fin.readline()[:-1]  # omit newline character
		all_data.append([0.0, -1, float(orig_score), 1.0000, orig_seq])
		# Get info for rest of sequences
		for line in fin:
			if (line.startswith(">")):
				seq_info = line[1:-1].split(", ")
				all_data.append([
								float(seq_info[0].split("=")[-1]),
								int(seq_info[1].split("=")[-1]),
								float(seq_info[2].split("=")[-1]),
								float(seq_info[3].split("=")[-1])
								])
				all_data[-1].append(fin.readline()[:-1])
	return pd.DataFrame(all_data, columns=["T", "sample", "score", "seq_recovery", "sequence"])

def add_name_list_to_MPNN_df(df, namebase, original_name):
	"""
	Generates a sequential list of names based on a namebase string + '_rank#' then sets the index to those names. Intended to be applied to already-sorted DFs (by score) to store the ordering in the FASTA output name. Assumes the starting sequence is the lowest rank (highest score) and labels it as original_name to distinguish it from MPNN outputs.
	"""
	namelist = [f"{namebase}_rank{i}" for i in range(1,len(df))]
	namelist.append(original_name)

	if (df.iloc[-1]["sample"] != -1):
		print("Original sequence is not highest MPNN score! Cannot order. Exiting...")
		return df

	# Add name list to column and change index to names (match regular fasta df)
	df["name"] = namelist
	return df.set_index("name")

def fasta_to_df(fname):
	"""
	Reads in a standard fasta file and returns as pandas df with columns [name, sequence]
	"""
	all_data = []
	with open(fname, "r") as fin:
		# Read each name/sequence pair and append to list. Skip comments & blank lines
		for line in fin:
			if (line.startswith("#")):
				continue
			if (line.startswith(">")):
				seq = fin.readline()
				all_data.append([line[1:-1], seq[:-1]])
			elif not line:
				continue
	return pd.DataFrame(all_data, columns=["name", "sequence"]).set_index("name")

def write_df_to_fasta(df, fname):
	"""
	Write sequence DF to standard fasta file. Supply multi-row DFs to get combined FASTA files or single DF lines to do individual outputs.
	"""
	with open(fname, "w") as fout:
		for name in df.index:
			fout.write(f">{df.loc[name].name}\n{df.loc[name].sequence}\n")
	print(f"Successfully wrote fasta file {fname} from dataframe.")

def add_loop(seq, loop_seq):
	"""
	Replaces the '/' in the MPNN sequence with the provided loop_seq string.
	"""
	return seq.replace("/", loop_seq)

def calc_num_residue_df(df, res_1letter):
	"""
	Calculate number of a given residue type (1 letter code) for all sequences in a dataframe and append as a new column called num_{res_1letter}.
	"""
	df[f"num_{res_1letter}"] = [len(seq.split(res_1letter))-1 for seq in df["sequence"]]

def calc_seq_identity_array(seq1, seq2):
	"""
	Takes two sequences of equal length and returns a list of ints corresponding to their per-position identity.
	"""
	if (len(seq1) != len(seq2)):
		print("The two sequences must be of equal length!! Returning False...")
		return False
	return [1 if seq1[i] == seq2[i] else 0 for i in range(len(seq1))]
