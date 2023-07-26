from hatch import *
from rand_walk import *
from greedy import *


# --------------------------------------------------------------------------- #
# Run clique finding algorithm
# --------------------------------------------------------------------------- #

def run_eclip_seed_cliques():
	global seeds, fnames
	seeds, fnames = [], []
	
	for eclip_id in get_eCLIP_ids():
		m = read_eclip_file(f"{LD}eCLIP/{eclip_id}.bed.gz")
		
		# get bin ids with most peaks in them
		top_bins = list(reversed(np.argsort(m)))[0:5]
		
		seeds.append(top_bins)
		fnames.append(eclip_id)
		
	print(f"\nread eCLIP data total = {len(fnames)}")
	
	results = mp.Pool(10).map(wroker_eclip, [i for i in range(len(fnames))])
	

def wroker_eclip(i):
	pi = do_random_walk(a, alpha, seeds[i])
	write_pi(f"{fout}/{fnames[i]}_{alpha}.txt", pi)


def read_eclip_file(fname):
	m = [0]*nbins
	with gzip.open(f"{fname}",'rb') as file_in:
		for line in file_in:
			words = line.decode('utf8').rstrip().split("\t")
			chr_name = words[0]

			if chr_name not in chr_map:
				continue

			bin = get_bin_from_chr_loc(chr_name, int(words[1]))
			m[bin] += 1
		
	return m


# --------------------------------------------------------------------------- #
# Run clique finding algorithm
# --------------------------------------------------------------------------- #

def run_random_seeds():
	results = mp.Pool(10).map(wroker_random, [i for i in range(100)])

def wroker_random(i):
	rand_bins = list(np.random.choice(nbins, 5, replace=False))
	pi = do_random_walk(a, alpha, rand_bins)
	write_pi(f"{fout}/{i}_{alpha}.txt", pi)


def find_random_seed():
	global e_id, top_bins_weight

	for eclip_id in get_eCLIP_ids():
		m = read_eclip_file(f"{LD}eCLIP/{eclip_id}.bed.gz")
		
		top_bins = list(reversed(np.argsort(m)))[0:5]
		
		top_bins_weight = find_clique_weight(top_bins, a)
		e_id = eclip_id


		results = mp.Pool(10).map(worker_random_enhanced, [i for i in range(20)])
	

def worker_random_enhanced(i):
	draws, dw = [], []
	# draw until you find a close enough seed
	for tt in range(50000):
		rand_bins = list(np.random.choice(nbins, 5, replace=False))
		rand_weight = find_clique_weight(rand_bins, a)
		
		
		if rand_weight/top_bins_weight >= 1:
			draws.append(rand_bins); dw.append(rand_weight)
			
	rand_bins = draws[list(np.argsort(dw))[0]]
	rand_weight = dw[list(np.argsort(dw))[0]]
	
	pi = do_random_walk(a, alpha, rand_bins)
	write_pi(f"{fout}/{e_id}_{alpha}_{i}.txt", pi)

	return 0
		

if __name__ == "__main__":
	pass

