from base import *

# --------------------------------------------------------------------------- #
# Helper functions
# --------------------------------------------------------------------------- #

def get_clique_total_weight(a, order):
	return sum([a[i][j] for (i,j) in it.combinations(order,2)])


def read_greedy_order(fname):
	return [(int(line.rstrip()),1) for line in open(fname)]


# write down the orders we found
def write_down_top_orders(a, top_orders, fname):
	for i in range(len(top_orders)):
		fo = open(fname, "w")
		for b in top_orders[i][0]:
			fo.write(f"{b}\n")
		fo.close()


# initialize the seed clique
def seed_clique_at(a, seed, used_bins):
	print(f"len(seed) = {len(seed)}")
	
	top_orders = [[seed[:], get_clique_total_weight(a, seed), "seeder"]]
	
	return top_orders
	
# call this function form the outside
def entry_grow_clique_greedily(a, seed, fname):
	start_time = time.time()
	
	used_bins = find_zero_and_non_zero_rows(a)[0]
	
	top_orders = grow_clique_greedily(a, seed, used_bins)
	
	print_time_taken(start_time, f"total running greedy = ")
	
	write_down_top_orders(a, top_orders, f"{fname}_gr")
	

# --------------------------------------------------------------------------- #
# Main recursive function
# --------------------------------------------------------------------------- #

def grow_clique_greedily(a, seed, used_bins):
	ZZ_MAX_CLIQUE_SIZE = 100
	
	top_orders = seed_clique_at(a, seed, used_bins)
	
	# 2. Expand the clique greedily
	while True:
		start_time = time.time()
		
		if len(top_orders[0][0]) == ZZ_MAX_CLIQUE_SIZE:
			print(f"max size reached\n")
			return top_orders
		
		possible_new_orders = []
		
		# for every good order we have so far, try expanding it
		for o in top_orders:
			order = o[0]
			delta = []
			
			# check every remaining bin to possibly add to the current clique
			for k in range(len(a)):
				if k in order or k in used_bins:
					continue
					
				scores = [a[i][k] for i in order]
							
				delta.append([k, score, sum(scores)])
				
			# sort them all to find the best increases for the given order
			delta.sort(key=lambda x: float(x[1]), reverse=True)
			
			# select the top orders for future expansion
			N = 1; M = 0
			# the top N + randomly drawn M from (N,7N)
			selected = list(range(N)) + random.sample(list(range(N,7*N)), M)
			for i in selected:
				new_order = order.copy()
				new_order.append(delta[i][0])
				new_score = o[1] + delta[i][2]
				
				possible_new_orders.append([new_order, new_score, "parent_id"])

		# sort the possible new orders by their score
		possible_new_orders.sort(key=lambda x: float(x[1]), reverse=True)
		
		# select the top orders for future expansion
		X1 = 1; X2 = 0; selected = list(range(len(possible_new_orders)))
		if len(possible_new_orders) > 100:
			selected = list(range(X1)) + random.sample(list(range(X1,len(possible_new_orders))), X2)

		top_orders = [possible_new_orders[i] for i in selected]
		
		print_time_taken(start_time, f"num_bins = {len(top_orders[0][0])}")
		

if __name__ == "__main__":
	pass
