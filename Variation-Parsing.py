#!/usr/bin/env python3

def slidingWindow(sequence,winSize,step):
	"""Returns a generator that will iterate through
	the defined chunks of input sequence.  Input sequence
	must be iterable.
	From scipher.wordpress.com"""
	# Verify the inputs
	try: it = iter(sequence)
	except TypeError:
		raise Exception("**ERROR** sequence must be iterable.")
	if not ((type(winSize) == type(0)) and (type(step) == type(0))):
		raise Exception("**ERROR** type(winSize) and type(step) must be int.")
	if step > winSize:
		raise Exception("**ERROR** step must not be larger than winSize.")
	if winSize > len(sequence):
		raise Exception("**ERROR** winSize must not be larger than sequence length.")

	# Pre-compute number of chunks to emit
	numOfChunks = int(((len(sequence)-winSize)/step)+1)
	# Do the work
	return numOfChunks
		
def make_whisker_plot(comparison, outname, style):
	#Displays a boxplot of the variance data generated for the variance
	#DNA compisition elements
	import seaborn as sns
	import matplotlib.pyplot as plt
	import pandas as pd
	#Due to the way that seaborn reads in the dataframe. The x-axis variable
	#must be repeated for each entry. The the y-axis data presented in a single
	#column
	#Three different blocks of text to address the different name variables provided
	#in the dataframes. Also used to rename the y-axis variables.
	if style == 'variance':
		df = pd.DataFrame(columns=['name','variance'])
		index_value = 0
		#Iterates through each set of values and adds a single row to dataframe
		#in order for correct seaborn format
		for i in comparison:
			for x in comparison[i]:
				df.loc[index_value] = [i, x]
				index_value += 1
		ax1 = sns.boxplot(x="name", y="variance", data=df)
		#Option for including points within figure - not recommended for TETRA
		#scores
		#Adds individual data points, transpartent alpha=0.2
		#Grey in color
		ax1 = sns.stripplot(x="name", y="variance", data=df, jitter=True, alpha=.2, color="grey")
		plt.xticks(rotation=90)
		fig1 = ax1.get_figure()
		#Set Y-axis size
		ax1.set_ylim([0,0.05])
		fig1.set_size_inches(10, 10)
		fig1.savefig(str(outname))
	if style == 'length':
		df = pd.DataFrame(columns=['name','length'])
		index_value = 0
		for i in comparison:
			for x in comparison[i]:
				df.loc[index_value] = [i, x]
				index_value += 1
		ax2 = sns.boxplot(x="name", y="length", data=df)
		ax2 = sns.stripplot(x="name", y="length", data=df, jitter=True, alpha=.2, color="grey")
		plt.xticks(rotation=90)
		#Set Y-axis size
		fig2 = ax2.get_figure()
		ax2.set_ylim([0,0.1])
		fig2.set_size_inches(10, 10)
		fig2.savefig(str(outname))
	if style == 'identity':
		df = pd.DataFrame(columns=['name','identity'])
		index_value = 0
		for i in comparison:
			for x in comparison[i]:
				df.loc[index_value] = [i, x]
				index_value += 1
		ax3 = sns.boxplot(x="name", y="identity", data=df)
		ax3.set_ylim([75,100])
		ax3 = sns.stripplot(x="name", y="identity", data=df, jitter=True, alpha=.2, color="grey")
		plt.xticks(rotation=90)
		fig3 = ax3.get_figure()
		fig3.set_size_inches(10, 10)
		fig3.savefig(str(outname))

def tet_clean(s):
    """ 
    Original code from Leighton Pritchard, leighton.pritchard@hutton.ac.uk
    redistributed and modified it under the terms of the GNU General 
    Public License as published by the Free Software Foundation, either 
    version 3 of the License, or (at your option) any later version.
    Checks that a passed string contains only unambiguous IUPAC nucleotide
    symbols. We are assuming that a low frequency of IUPAC ambiguity symbols
    doesn't affect our calculation.
    """
    if not len(set(s) - set('ACGT')):
        return True
    return False


def calc_org_tetra(fn, org):
    """ 
    Original code from Leighton Pritchard, leighton.pritchard@hutton.ac.uk
    redistributed and modified it under the terms of the GNU General 
    Public License as published by the Free Software Foundation, either 
    version 3 of the License, or (at your option) any later version.
    Calculate the tetranucleotide frequencies
    for each sequence, on each strand, and follow Teeling et al. (2004)
    in calculating a corresponding Z-score for each observed
    tetranucleotide frequency, dependent on the mono-, di- and tri-
    nucleotide frequencies for that input sequence.
    """
    import collections
    from Bio import SeqIO
    import math
    org_tetraz = {}
    # For the Teeling et al. method, the Z-scores require us to count
    # mono, di, tri and tetranucleotide sequences
    monocnt, dicnt, tricnt, tetracnt = (collections.defaultdict(int),
                                        collections.defaultdict(int),
                                        collections.defaultdict(int),
                                        collections.defaultdict(int))

    for s in [str(fn).upper(),
              str(fn.reverse_complement()).upper()]:
        # Since the Teeling et al. algorithm requires us to consider
        # both strand orientations, monocounts are easy
        monocnt['G'] += s.count('G')
        monocnt['C'] += s.count('C')
        monocnt['T'] += s.count('T')
        monocnt['A'] += s.count('A')
        # For di, tri and tetranucleotide counts, we loop over the
        # sequence and its reverse complement, until we're near the end:
        for i in range(len(s[:-4])):
            di, tri, tetra = s[i:i+2], s[i:i+3], s[i:i+4]
            dicnt[str(di)] += 1
            tricnt[str(tri)] += 1
            tetracnt[str(tetra)] += 1
        # We clean up the straggling bit at the end:
        tricnt[str(s[-4:-1])] += 1
        tricnt[str(s[-3:])] += 1
        dicnt[str(s[-4:-2])] += 1
        dicnt[str(s[-3:-1])] += 1
        dicnt[str(s[-2:])] += 1
    # Following Teeling (2004), we calculate expected frequencies for each
    # tetranucleotide; we ignore ambiguity symbols
    tetra_exp = {}
    for t in [tet for tet in tetracnt if tet_clean(tet)]:
        tetra_exp[t] = 1.*tricnt[t[:3]]*tricnt[t[1:]]/dicnt[t[1:3]]
    # Following Teeling (2004) we approximate the std dev of each
    # tetranucleotide
    tetra_sd = {}
    for t, exp in tetra_exp.items():
        den = dicnt[t[1:3]]
        tetra_sd[t] = math.sqrt(exp * (den - tricnt[t[:3]]) * \
                                    (den - tricnt[t[1:]]) / (den * den))
    # Following Teeling (2004) we calculate the Z-score for each
    # tetranucleotide
    tetra_z = {}
    for t, exp in tetra_exp.items():
        try:
            tetra_z[t] = (tetracnt[t] - exp)/tetra_sd[t]
        except ZeroDivisionError:
            # We hit a zero in the estimation of variance
            zeroes = [k for k,v in tetra_sd.items() if v == 0]
            tetra_z[t] = 1 / (dicnt[t[1:3]] * dicnt[t[1:3]])
    org_tetraz[org] = tetra_z
    return org_tetraz

def calc_tetra_corr(tetra_z):
    """ 
    Original code from Leighton Pritchard, leighton.pritchard@hutton.ac.uk
    redistributed and modified it under the terms of the GNU General 
    Public License as published by the Free Software Foundation, either 
    version 3 of the License, or (at your option) any later version.
    Calculate Pearson correlation coefficient from Z scores for each
    tetranucleotide. 
    """
    from scipy.stats.stats import pearsonr
    import numpy as np
    corrs = []
    orgs = sorted(tetra_z.keys())
    #Iterates through the dictionary in a unidirectional manner to
    #generate pairwise comaprison
    for idx, o1 in enumerate(orgs[:-1]):
        for o2 in orgs[idx+1:]:
    #Modified from the original script to use numpy arrays and
    #pearson correlation module
            z1_list = []
            z2_list = []
            for k in sorted(tetra_z[o1].keys()):
                try:
                    z2_list.append(tetra_z[o2][k])
                    z1_list.append(tetra_z[o1][k])
                except KeyError:
                    continue
            z1_array = np.asarray(z1_list)
            z2_array = np.asarray(z2_list)
            corrs.append(pearsonr(z1_array,z2_array)[0])
    return corrs

def run_stats(names, comparison):
	import pandas as pd
	from scipy import stats
	df = pd.DataFrame(columns=['name','variance'])
	index_value = 0
	#Iterates through each set of values and adds a single row to dataframe
	#in order for correct seaborn format
	for i in comparison:
		for x in comparison[i]:
			df.loc[index_value] = [i, x]
			index_value += 1
	#After clean up, use Wilks-Shapiro test to determine data normalcy
	normal_data = []
	for n in names:
		print(n, stats.shapiro(df['variance'][df['name'] == n]))
		#If p-value greater than alpha 0.05 accept null hypothesis
		#Data set is normally distributed
		if float(stats.shapiro(df['variance'][df['name'] == n])[1]) > 0.05:
			normal_data.append(n)
	print(str(normal_data))
	#For normal data sets determine the T-test significance test
	if len(normal_data) > 1:
		for i,x in enumerate(normal_data):
			for j,y in enumerate(normal_data):
				if i > j:
					print(str(x), str(y), stats.ttest_ind(df['variance'][df['name'] == x], df['variance'][df['name'] == y], equal_var=False))
	#For data sets that do not adhere to normalcy, determine Mann-Whitney Ranked Sums test	
	for i,x in enumerate(names):
		for j,y in enumerate(names):
			if i > j:
				print(str(x), str(y), stats.ranksums(df['variance'][df['name'] == x], df['variance'][df['name'] == y]))
	return


def main():
	from Bio import SeqIO
	import argparse
	import numpy as np
	import pandas as pd
	from Bio.SeqUtils import GC
	import subprocess
	import os.path
	from scipy import stats
	parser = argparse.ArgumentParser(description="Calculate variance of a set of genomes \
		for GC content, tetranucleotide variation, internal redundancy, and global coverage \
		value consistency")
	parser.add_argument('-i', '--input', help="Comma separate list of genome set filenames. \
		Must include at least 1 set of genomes. Recommend calculcation on >> 1 genome/ \
		If multiple sets provide a box and whisker output figure will be calculated")
	#Genome set format name of genomes before extension. Genome name will need to be same
	#between input file types
	parser.add_argument('-f', '--fastaext', help="File extension used for FASTA files")
	parser.add_argument('--runGC',action='store_true', help="Compute GC variance")
	parser.add_argument('--runTETRA',action='store_true', help="Compute tetranucleotide variance")
	parser.add_argument('--runRedundancy',action='store_true', help="Compute NUCmer redundancy variance")

	args = parser.parse_args()
	arg_dict = vars(args)

	if arg_dict['runGC'] == True:
		gc_comparison = {}
		#Populate dictionary that will be used for variance analysis
		#Splits genome set file names
		for i in arg_dict['input'].split(","):	
			#Each line represents a genome name
			if os.path.exists(str(i)+".gcvar") == False:
				outgcfile = open(str(i)+".gcvar", "w")
				for line in open(i, "r"):
					line = line.strip()
					genome_gc = []
					#Parse through each contig. Splitting each contig into
					#2000bp fragments with a 500bp step
					for record in SeqIO.parse(open(str(line)+"."+str(arg_dict['fastaext']), "r"), "fasta"):
						if len(record.seq) < 2000:
							continue
						numOfChunks = slidingWindow(record.seq,2000,500)
						for t in range(0,numOfChunks*500,500):
							genome_gc.append(GC(record.seq[t:t+500]))
					#Create a numpy array and calculate variance for all contig chunks
					#For entire genome
					gc_var = np.var(pd.DataFrame(genome_gc), axis=0)
					#Numpy array is a (1,2) array. Grab the second value as the
					#first value is 0 -- unclear why
					for x in np.nditer(gc_var):
						if x > 0:
							try:
								gc_comparison[i].append(float(x))
								outgcfile.write(i+"\t"+str(x)+"\n")
							except KeyError:
								gc_comparison[i] = [float(x)]
								outgcfile.write(i+"\t"+str(x)+"\n")
			#If the data has been calculated previously and stored in a checkpoint file
			#Process data to be run through stats and plot making
			else:
				for line in open(str(i)+".gcvar", "r"):
					line = line.strip()
					data = line.split("\t")
					try:
						gc_comparison[data[0]].append(float(data[1]))
					except KeyError:
						gc_comparison[data[0]] = [float(data[1])]

		names = arg_dict['input'].split(",")
		run_stats(names, gc_comparison)
		make_whisker_plot(gc_comparison, "gc_variance.svg", "variance")

	if arg_dict['runRedundancy'] == True:
		#Performs analysis using NUCmer to identify repeat regions within a genome
		red_len_comparison = {}
		red_ident_comparison = {}
		for i in arg_dict['input'].split(","):
			for line in open(i, "r"):
				line = line.strip()
				total_len = 0
				redundant_lengths = []
				redundant_identity = []
				#Determine the total length of each genome
				for record in SeqIO.parse(open(str(line)+"."+str(arg_dict['fastaext']), "r"), "fasta"):
					total_len += len(record.seq)		
				#Create NUCmer output file name
				coords_out = str(line)+".coords"
				previously_matched = []
				#Only run NUCmer is .coords output file is absent
				if os.path.exists(coords_out) == False:
					fasta_in = str(line)+"."+str(arg_dict['fastaext'])
					subprocess.call(["nucmer", "--maxmatch", "-p", line, fasta_in, fasta_in])
					delta_in = str(line)+".delta"
					#Open .coords outfile
					outfile = open(coords_out, "w")
					#Store .coords data in file
					subprocess.call(["show-coords", "-T", "-H", delta_in],stdout=outfile)
					outfile.close()
				#Parse the coordinate data in each .coords file for each genome
				for coords in open(coords_out, "r"):
					coords = coords.strip()
					data = coords.split("\t")
					#Create a unique string for each query and ref representing regions of overlap
					query = str(data[7])+"-"+str(min(data[0],data[1]))+"-"+str(max(data[0],data[1]))+"-"+data[6]
					ref = str(data[8])+"-"+str(min(data[2],data[3]))+"-"+str(max(data[2],data[3]))+"-"+data[6]
					#If the identified region is a self-hit, the same region identified as query and reference
					#skip
					if data[7] == data[8] and data[0] == data[2]:
						continue
					#Calculate all lengths where %identity is >=97
					if query not in previously_matched and float(data[6]) >=97:
						redundant_lengths.append(float(min([data[4],data[5]]))/float(total_len))
						redundant_identity.append(float(data[6]))
						previously_matched.append(query)
						previously_matched.append(ref)
				#Check to see if any regions in genome were in repeat regions before proceeding
				if len(redundant_lengths) > 0 and len(redundant_identity) > 0:
					try:
						#Exclude genomes with excessively high >50% of length in repeat regions
						#This can occur when there are multiple copies (3+) of the same repeat
						#Each instance counts separately 
						if sum(redundant_lengths) < 0.5:
							red_len_comparison[i].append(sum(redundant_lengths))
						red_ident_comparison[i].append(sum(redundant_identity)/len(redundant_identity))
					except KeyError:
						red_len_comparison[i] = [sum(redundant_lengths)]
						red_ident_comparison[i] = [sum(redundant_identity)/len(redundant_identity)]
				
		names = arg_dict['input'].split(",")
		run_stats(names, red_len_comparison)
		make_whisker_plot(red_len_comparison, "redundancy_length.svg", "length")
		make_whisker_plot(red_ident_comparison, "redundant_avg_identity.svg", "identity")

	if arg_dict['runTETRA'] == True:
		tet_comparison = {}
		for i in arg_dict['input'].split(","):
			if os.path.exists(str(i)+".tetravar") == False:
				outtetrafile = open(str(i)+".tetravar", "w")
				for line in open(i, "r"):
					line = line.strip()
					chunk_zscore = {}
					#Parse through each contig. Splitting each contig into
					#10000bp fragments with a 5000bp step
					#To change values, modify below
					for record in SeqIO.parse(open(str(line)+"."+str(arg_dict['fastaext']), "r"), "fasta"):
						#Here
						if len(record.seq) < 10000:
							continue
						#Here
						numOfChunks = slidingWindow(record.seq,10000,5000)
						#Here
						for t in range(0,numOfChunks*5000,5000):
							#Here
							chunk_zscore.update(calc_org_tetra(record.seq[t:t+5000],t))
					tet_var = np.var(pd.DataFrame(calc_tetra_corr(chunk_zscore)), axis=0)
					try:
						for x in np.nditer(tet_var):
							if x > 0:
								try:
									tet_comparison[i].append(float(x))
									outtetrafile.write(i+"\t"+str(x)+"\n")
								except KeyError:
									tet_comparison[i] = [float(x)]
									outtetrafile.write(i+"\t"+str(x)+"\n")
					except ValueError:
						continue
			#If the data has been calculated previously and stored in a checkpoint file
			#Process data to be run through stats and plot making
			else:
				for line in open(str(i)+".tetravar", "r"):
					line = line.strip()
					data = line.split("\t")
					try:
						tet_comparison[data[0]].append(float(data[1]))
					except KeyError:
						tet_comparison[data[0]] = [float(data[1])]
				
		names = arg_dict['input'].split(",")
		run_stats(names, tet_comparison)
		make_whisker_plot(tet_comparison, "tetra_variance.svg", "variance")


if __name__ == "__main__":
	main()