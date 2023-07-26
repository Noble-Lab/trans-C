from base import *


# --------------------------------------------------------------------------- #
# Convert from hic to mcool
# --------------------------------------------------------------------------- #

def convert_hic_to_mcool(idd):
	cmd = f"hicConvertFormat -m {idd}.hic -o {idd} --inputFormat hic --outputFormat cool"



# --------------------------------------------------------------------------- #
# Convert from mcool to npy
# --------------------------------------------------------------------------- #

def convert_mcool_to_npy(fmcool, fname_out, res):
	import cooler
	
	start_time = time.time()
	
	c = cooler.Cooler(f"{fmcool}::/resolutions/{res}")
	
	x = c.matrix(balance=False)
	
	nbins = len(x)
	print(f"nbins = {nbins}")
	
	# this saves TONS of time
	z = x[0:nbins]

	a = np.zeros((nbins, nbins))
	
	for i in range(nbins):
		y = z[i]

		for j in range(i,nbins):
			a[i][j] = y[j]
			a[j][i] = a[i][j]
			
	
	print_time_taken(start_time, f"all for res = {res}")
	np.save(fname_out, a)
	return a
	

# --------------------------------------------------------------------------- #
# Convert intact Hi-C from hic to npy
# --------------------------------------------------------------------------- #

def convert_intact_hic_to_npy(fname, fname_out, res):
	a = np.zeros((nbins, nbins))
	
	for chr1 in chr_map:
		for chr2 in chr_map:
			# use hic-straw to extract the counts per pair of chromosomes
			result = hicstraw.straw('observed', 'NONE', f"{fname}", chr1, chr2, 'BP', res)
			
			for i in range(len(results)):
				# hic-straw returns the beggining of the interval
				bin1 = get_bin_from_chr_loc(chr1, result[i].binX + 1)
				bin2 = get_bin_from_chr_loc(chr2, result[i].binY + 1)
				
				a[bin1][bin2] = result[i].counts
				
				# the matrix is symmetric
				a[bin2][bin1] = a[bin1][bin1]
	
	np.save(fname_out, a)	
	return a
	


# --------------------------------------------------------------------------- #
# Convert from .pairs.gz to npy
# --------------------------------------------------------------------------- #

def convert_pairs_to_npy(fname, fname_out):
	m = np.zeros((nbins, nbins))
	
	with gzip.open(f"{fname}.pairs.gz",'rt') as f:
		for line in f:
			words = line.rstrip()
						
			# skip comments
			if words[0] == "#":
				continue
			
			num_lines += 1
			
			words   = line.split("\t")
			read_id = words[0]
			chr1    = words[1]
			pos1	= int(words[2])
			chr2    = words[3]
			pos2	= int(words[4])
			
			# skip non chromosomes, aka special contigs
			if chr1 not in chr_map or chr2 not in chr_map:
				continue
			
			bin1 = get_bin_from_chr_loc(chr1, pos1)
			bin2 = get_bin_from_chr_loc(chr2, pos2)
			
			m[bin1][bin2] += 1
			
			# The matrix is symmetric !!!
			m[bin2][bin1] += 1
	
	print(f"read {fname}.pairs.gz with num_lines = {num_lines}")
	np.save(fname_out, m)
	return m





# --------------------------------------------------------------------------- #
# Ice Normalize a matrix
# --------------------------------------------------------------------------- #

def ice_normalize_matrix(a, fname, bias_name=""):
	from iced import normalization
	from iced import filter
	
	start_time = time.time()
	counts = filter.filter_low_counts(a, percentage=0.05)
	(normed, bias) = normalization.ICE_normalization(counts, output_bias=True)
	
	# write the biases down
	if len(bias_name) > 1:
		fo = open(bias_name, "w")
		for i in bias:
			fo.write(str(i[0]) + "\n")
	
	# write the ICED matrix down
	#write_matrix(normed, fname, False)
	np.save(fname, normed)
	print_time_taken(start_time, f"icing matrix {fname}")
	return normed


# --------------------------------------------------------------------------- #
# Rescales the heaviest weights in a matrix above given threshold
# --------------------------------------------------------------------------- #

def zero_cis_clip_extreme(a, fname_out):
	# zero all cis contatcs
	for i in range(nbins):
		for j in range(i, nbins):
			# zero the cis-contacts
			if bin_map[i][0] ==  bin_map[j][0]: 
				a[i][j] = 0
				a[j][i] = 0
			
	# find the threshold to clip
	vals = a.flatten()
	sorted_vals = sorted(vals)
	thresh = sorted_vals[int(0.99999*len(vals))]
	
	# clip the trans-contatcs above the threshold
	for i in range(nbins):
		for j in range(i, nbins):
			# log scale trans-contcts that are greater than the threshold value
			if a[i][j] > thresh:
				a[i][j] = thresh + math.log(1 + a[i][j] - thresh)
				a[j][i] = a[i][j]
	
	np.save(fname_out, a)
	return a



# --------------------------------------------------------------------------- #
# Create bin maps
# --------------------------------------------------------------------------- #

# Creates a bin map. uses chr_sizes. assumes chr are in order. keeps chr_id in the first column
def create_bin_map(fpath):
	# carefully select the bin size and the chromosome lengths and set the name of the organism
	(chr_length, chr_names) = read_chr_lengths()
	
	for bin_size in [10000, 50000, 100000, 250000, 500000, 1000000]:

		bin, start, fo = 0, 0, open(f"{fpath}bin_map_{bin_size}.bed", "w")
	
		for chr in chr_names:
			while True:
				end = start + bin_size
			
				fo.write(f"{chr}\t{start}\t{end}\t{bin}\n")
			
				if end < chr_length[chr]:
					start += bin_size
					bin   += 1
				else:
					start = 0
					bin   += 1
					break


if __name__ == "__main__":
	pass


