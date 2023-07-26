from hatch import *
from rand_walk import *
from greedy import *

	

def TF_matched_random_seed():
	global tf_ids
	tf_ids = get_tf_ids("GM12878")
		
	results = mp.Pool(10).map(TF_worker_matched_random, [i for i in range(len(tf_ids))])


def TF_worker_matched_random(i):
	tf_id = tf_ids[i]
	m = read_tf_bed(tf_id)
	
	top_bins = list(reversed(np.argsort(m)))[0:5]
	top_bins_weight = find_clique_weight(top_bins, a)
	
	for num_samples in range(10):
		# draw until you find a close enough seed
		draws, dw = [], []
		for tt in range(50000):
			rand_bins = list(np.random.choice(nbins, 5, replace=False))
			rand_weight = find_clique_weight(rand_bins, a)
		
			if rand_weight/top_bins_weight >= 1:
				draws.append(rand_bins); dw.append(rand_weight)
		
		if len(dw) > 0:
			rand_bins = draws[list(np.argsort(dw))[0]]
			rand_weight = dw[list(np.argsort(dw))[0]]
		
			pi = do_random_walk(a, alpha, rand_bins)
			write_pi(f"{fout}/{tf_id}_{alpha}_{num_samples}.txt", pi)
		

def TF_purely_random_seed():
	results = mp.Pool(10).map(TF_worker_purely_random, [i for i in range(50)])
	

def TF_worker_purely_random(i):
	rand_bins = list(np.random.choice(nbins, 5, replace=False))
	pi = do_random_walk(a, alpha, rand_bins)
	write_pi(f"{fout}/{i}_{alpha}_mrw.txt", pi)
	
	
def TF_seed():
	global tf_ids
	tf_ids = get_tf_ids("GM12878")
		
	results = mp.Pool(18).map(TF_worker_seed, [i for i in range(len(tf_ids))])
	

def TF_worker_seed(i):
	tf_id = tf_ids[i]
	m = read_tf_bed(tf_id)
	
	top_bins = list(reversed(np.argsort(m)))[0:5]
	pi = do_random_walk(a, alpha, top_bins)
	write_pi(f"{fout}/{tf_id}_{alpha}.txt", pi)


def run_tf_seed_cliques():
	global seeds, fnames
	seeds, fnames = [], []
	
	for tf_id in get_tf_ids("GM12878"):
		m = read_tf_bed(tf_id)

		# get bin ids with most peaks in them
		top_bins = list(reversed(np.argsort(m)))[0:5]
		
		seeds.append(top_bins)
		fnames.append(tf_id)

	print(f"\nread all tf_ids total = {len(fnames)}")
	
	results = mp.Pool(10).map(wroker_tf1, [i for i in range(len(fnames))])




if __name__ == "__main__":
	pass

