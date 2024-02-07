from base import *

import pyBigWig

# --------------------------------------------------------------------------- #
# TSA-seq data
# --------------------------------------------------------------------------- #

# read the TSA-seq data
def read_tsa_seq(fname):
	bw = pyBigWig.open(f"{D}tsa_seq/{fname}.bw")

	tsa = [0]*nbins
	for b in range(nbins):
		(chr, start, end) = bin_map[b]
		
		# extract the value into our array
		tsa[b] = float(bw.stats(chr, start, end)[0])
		
		# NOTE: this doesn't handle NANs! will return NAN
		# tsa[b] = np.average(bw.values(f"{chr}", start, end))

	return tsa

# read the cliques pvalues
def read_pvals(fname):
	hsh = {}
	for line in open(f"{D}pvals/{fname}_pvals.txt"):
		words = line.rstrip().split("\t")
		hsh[words[0]] = float(words[3])
	
	return hsh

def read_tsas(type):
	hsh = {}
	for line in open(f"{T}tsa_{type}.txt"):
		words = line.rstrip().split("\t")
		hsh[words[2]] = float(words[0])

	return hsh



# TF cliques analysis
def run_tf_sorted(fname):
	NUM_BINS = 40
	tsa = read_tsa_seq(fname)
	
	fo = open(f"{fout}tsa.txt", "w")
	
	for tf_id in get_tf_ids("GM12878"):
		m = read_tf_bed(tf_id)

		# get bin ids with most peaks in them
		sorted_bins = list(reversed(np.argsort(m)))
		top_bins = sorted_bins[0:NUM_BINS]
		
		# find bins with fewest > 0 peaks
		for i in range(nbins):
			if sorted_bins[i] == 0:
				min_val - sorted_bins[i-1]
				min_val_bins = [j for j in range(i) if sorted_bins[j] == min_val]
				bottom_bins = min_val_bins[-NUM_BINS:]
				break
		
		
		top_tsa, bottom_tsa = [], []
		for i in range(len(top_bins)):
			top_tsa.append(tsa[top_bins[i]])
			bottom_tsa.append(tsa[bottom_bins[i]])
			
		# write down for analysis
		print(f"np.average(top_tsa) = {np.nanmean(top_tsa)} np.average(bottom_tsa) = {np.nanmean(bottom_tsa)}\n")
		fo.write(f"{np.nanmean(top_tsa)}\t{np.nanmean(bottom_tsa)}\t{tf_id}\n")



def run_tf_clique(fname):
	hsh_pval = read_pvals("tf")
	hsh_names = read_tf_ids_to_gene_names()
	
	NUM_BINS = 40
	tsa = read_tsa_seq(fname)
	
	fo = open(f"{fout}tsa_tf.txt", "w")
	
	# for each TF
	for tf_id in get_tf_ids("GM12878"):
		m = read_tf_bed(tf_id)

		pi = read_pi(f"{P}Outs/TF_seed/{tf_id}.txt")
		clique = [x[0] for x in pi][0:NUM_BINS]
		
		arr_rank = [x[0] for x in pi]
		tf_rank = []
		for b in arr_rank:
			if m[b] > 0:
				tf_rank.append(b)
		
		# lowest ranked bins with > 0 peaks
		bottom = tf_rank[-NUM_BINS:]
		
		top_tsa, bottom_tsa = [], []
		for i in range(len(clique)):
			top_tsa.append(tsa[clique[i]])
			bottom_tsa.append(tsa[bottom[i]])
			
		print(f"np.average(top_tsa) = {np.nanmean(top_tsa)} np.average(bottom_tsa) = {np.nanmean(bottom_tsa)}")

		color = ("no" if hsh_pval[hsh_names[tf_id]] < 0.05 else "sig")
		fo.write(f"{np.nanmean(top_tsa)}\t{np.nanmean(bottom_tsa)}\t{hsh_names[tf_id]}\t{color}\n")



# eCLIP analysis
def run_eclip_sorted(fname):
	NUM_BINS = 40
	tsa = read_tsa_seq(fname)
	
	fo = open(f"{fout}tsa_eclip_sorted.txt", "w")
	
	for tf_id in get_eCLIP_ids():
		m = read_eclip_file(f"{LD}eCLIP/{tf_id}.bed.gz")
		
		# get bin ids with most peaks in them
		sorted_bins = list(reversed(np.argsort(m)))
		top_bins = sorted_bins[0:NUM_BINS]
		
		# find bins with fewest > 0 peaks
		for i in range(nbins):
			if sorted_bins[i] == 0:
				min_val - sorted_bins[i-1]
				min_val_bins = [j for j in range(i) if sorted_bins[j] == min_val]
				bottom_bins = min_val_bins[-NUM_BINS:]
				break
		
		top_tsa, bottom_tsa = [], []
		for i in range(len(top_bins)):
			top_tsa.append(tsa[top_bins[i]])
			bottom_tsa.append(tsa[bottom_bins[i]])
			
		fo.write(f"{np.nanmean(top_tsa)}\t{np.nanmean(bottom_tsa)}\t{tf_id}\n")


def run_eclip_clique():
	hsh_pval = read_pvals("eclip")
	hsh_names = read_eclip_ids_to_gene_names()
	
	NUM_BINS = 40
	tsa = read_tsa_seq()
	
	fo = open(f"{fout}tsa_eclip.txt", "w")
	
	for tf_id in get_eCLIP_ids():
		m = read_eclip_file(f"{LD}eCLIP/{tf_id}.bed.gz")

		pi = read_pi(f"{P}Outs/eCLIP_seed/{tf_id}.txt")
		clique = [x[0] for x in pi][0:NUM_BINS]
		
		arr_rank = [x[0] for x in pi]
		tf_rank = []
		for b in arr_rank:
			if m[b] > 0:
				tf_rank.append(b)
		
		bottom = tf_rank[-NUM_BINS:]
		
		top_tsa, bottom_tsa = [], []
		for i in range(len(clique)):
			top_tsa.append(tsa[clique[i]])
			bottom_tsa.append(tsa[bottom[i]])
			
		print(f"np.average(top_tsa) = {np.nanmean(top_tsa)} np.average(bottom_tsa) = {np.nanmean(bottom_tsa)}")
		
		color = ("no" if hsh_pval[hsh_names[tf_id]] < 0.05 else "sig")
		fo.write(f"{np.nanmean(top_tsa)}\t{np.nanmean(bottom_tsa)}\t{hsh_names[tf_id]}\t{color}\n")
		



# calculate p-value per file
def tsa_stats(fname):
	clique, bottom, count = [], [], 0
	for line in open(f"{T}{fname}.txt"):
		words = line.rstrip().split("\t")
		clique.append(float(words[0]))
		bottom.append(float(words[1]))
		if float(words[0]) > float(words[1]):
			count += 1
	
	print(f"count = {count} out of total = {len(clique)}")
	x = stats.mannwhitneyu(clique, bottom, alternative = "less")
	print(f"x.pvalue = {x.pvalue}")


def plot_rank():
	typ = "tf"
	hsh_pval = read_pvals(typ)
	hsh_tsas = read_tsas(typ)
	
	fo = open(f"{fout}ranks_{typ}.txt", "w")
	for tf_id in hsh_pval:
		if tf_id not in hsh_tsas:
			continue
		fo.write(f"{hsh_tsas[tf_id]}\t{hsh_pval[tf_id]}\t{tf_id}\n")
	


# entry for the analysis
def entry_start():
	global idd, organism, res, chr_length, chr_names, chr_map, bin_map, nbins
	organism = "human"
	res      = 100000
	
	(chr_length, chr_names, chr_map, bin_map, nbins) = read_fixed_inputs(organism, res)
	
	run_tf_clique()
	tsa_stats("tsa_tf")
	run_eclip_clique()
	tsa_stats("tsa_eclip")
	


# --------------------------------------------------------------------------- #
# Main Executable
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
	entry_start()


