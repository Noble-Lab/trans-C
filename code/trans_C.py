from base import *
from rand_walk import *

	
def read_user_inputs():
	print(f"Starting trans-C.\n")
	
	if len(sys.argv) < 8:
		print(f"insufficient arguments.\n")
	
	global fpath, a, res, seed, chrom_sizes, bin_map
	
	a = np.load(sys.argv[2])
	bin_map = create_bin_map(sys.argv[3])
	res = int(sys.argv[4])
	seed = read_seed(sys.argv[5])
	fpath = sys.argv[6])
	alpha = float(sys.argv[7])

	
	do_random_walk(a, gamma, start_locs)
	
	clean_up()
	
	print(f"\ntrans-C done.\n")

if __name__ == "__main__":
	read_user_inputs()

