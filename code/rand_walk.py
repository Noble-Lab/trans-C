from base import *

# --------------------------------------------------------------------------- #
# Do a Random Walk
# --------------------------------------------------------------------------- #	

# Do a random walk analytically
def do_random_walk(a, gamma, start_locs):
	start_time = time.time()
	nbins = len(a)
	
	# create the transitional matrix P
	p = np.zeros((nbins,nbins))
	
	# for every node i
	for i in range(nbins):
		outgoing = sum(a[i])
		
		# find the probability to move to node j
		for j in range(nbins):			
			# you can get via i -> j or you restart at j
			if outgoing > 0:
				p[i][j] = (1-gamma)*a[i][j]/outgoing + gamma*1.0*(j in start_locs)/len(start_locs)
			# special case: all zeros in that row i
			else:
				p[i][j] = gamma*1.0*(j in start_locs)/len(start_locs)

	# sanity test
	for i in range(nbins):
		break
		if not math.isclose(sum(p[i]),1):
			print(f"i = {i}; {sum(p[i])}")
			#sys.exit()
			pass

	
	# find the left eigenvector corresponding to eigenvalue 1
	w, vl = eig(p,left=True, right=False)
	
	# extract it (the real part)
	for i in range(len(w)):
		if math.isclose(w[i].real,1):
			pi = vl[:,i]
			break

	
	# normalize the vector
	pi = [float(i.real)/sum(pi).real for i in pi]
	
	
	st = print_time_taken(start_time, f"Doing RW via calculating eigenvectors for matrix of size {nbins}")
	
	fo = open(f"{T}analytical.txt", "a")
	fo.write(f"{st}\n")
	fo.close()
	
	return pi
	
	
def write_pi(fname, pi):
	fo = open(f"{fname}", "w")
	for i in reversed(np.argsort(pi)):
		fo.write(f"{i}\t{pi[i]}\n")
	fo.close()


if __name__ == "__main__":
	pass
	
	
