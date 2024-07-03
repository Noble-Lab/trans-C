from base import *
from tf import *

from scipy.spatial.distance import pdist,cdist,squareform

def process_proximity_matrix(fout):
	folder = f"{D}MERFISH/"
	save_data = fout

	# reads the data in
	experiment = []
	fid = open(folder+os.sep+r'genomic-scale.tsv','r')
	lines = np.array([ln[:-1].split('\t') for ln in fid if len(ln)>0])
	head = list(lines[0])
	experiment = np.concatenate([experiment,lines[1::2082,head.index('experiment number')].astype(int)])
	zxy = np.array(lines[1:,:3][:],dtype=float)
	dLAM = np.array(lines[1:,-1].astype(float))


	fid = open(folder+os.sep+r'genomic-scale-with transcription and nuclear bodies.tsv','r')
	lines = np.array([ln[:-1].split('\t') for ln in fid if len(ln)>0])
	head = list(lines[0])
	experiment = np.concatenate([experiment,lines[1::2082,head.index('experiment number')].astype(int)])
	dLAM = np.concatenate([dLAM,np.array(lines[1:,-3].astype(float))])
	zxy = np.concatenate([zxy,np.array(lines[1:,:3][:],dtype=float)])
	zxy = zxy.reshape([-1,2082,3])/1000 #transform to um
	dLAM = dLAM.reshape([-1,2082])/1000

	lens = [76, 80, 66, 63, 60, 55, 53, 48, 40, 43, 44, 44, 33, 30, 31, 30, 33, 33, 33, 33, 31, 31, 51]
	edges = [0]+list(np.cumsum(lens))
	ijs = []
	for i in range(len(lens)):
	    for j in range(len(lens)):
	        ijs.append((i,j))
	im_med = np.zeros([edges[-1],edges[-1]])
	cut_offs = [0.25,0.5,0.75,1]

	im_fr = np.zeros([edges[-1],edges[-1],len(cut_offs)])
	im_med_trans = []
	im_med_cis = []
	im_fr_trans = [[] for _ in cut_offs]
	im_fr_cis = [[] for _ in cut_offs]

	zxy_ = zxy
	for i,j in ijs:
	    arr = []
	    for st1 in [0,edges[-1]]:
	        for st2 in [0,edges[-1]]:
	            zxy1 = zxy_[:,st1+edges[i]:st1+edges[i+1]]
	            zxy2 = zxy_[:,st2+edges[j]:st2+edges[j+1]]
	            arr =arr+[cdist(zxy1[k],zxy2[k]) for k in range(len(zxy1))]
	    arr = np.array(arr)
	    im_med[edges[i]:edges[i+1],edges[j]:edges[j+1]]=np.nanmedian(arr,axis=0)
	    if i==j:
	        im_med_cis.append(np.nanmedian(arr[::2],axis=0))
	        im_med_trans.append(np.nanmedian(arr[1::2],axis=0))
	    for ic,cutoff in enumerate(cut_offs):
	        im_fr[edges[i]:edges[i+1],edges[j]:edges[j+1],ic] = 1.*np.sum(arr<cutoff,0)/np.sum(arr>-1,0)
	        if i==j:
	            im_fr_trans[ic].append(1.*np.sum(arr[1::2]<cutoff,0)/np.sum(arr[1::2]>-1,0))
	            im_fr_cis[ic].append(1.*np.sum(arr[::2]<cutoff,0)/np.sum(arr[::2]>-1,0))
	
	
	# save the matrix as a pickle object
	pickle.dump([im_med,im_fr,im_med_trans,im_med_cis,im_fr_trans,im_fr_cis,len(zxy)],
	        open(save_data+r'/mat_contact_IMR90_untreated.pkl','wb'))

# zero cis and normalzie trans
def get_proximity_weight(im_med, im_fr, clique_loci, arr_loci):
	weight_med, weight_fr, normlz = 0, 0, 0
	for i in clique_loci:
		for j in clique_loci:
			# skip cis-contacts
			if arr_loci[i][0] == arr_loci[j][0]:
				continue
			weight_med += im_med[i][j]
			weight_fr  += im_fr[i][j]
			normlz += 1

	return weight_fr/normlz


def get_merfish_loci_coordinates():
	hsh, arr, i  = {}, [], 0
	for line in open(f"{D}MERFISH/genomic-scale.tsv"):
		i += 1
		if i == 1:
			continue
		words = line.rstrip().split("\t")
		chr = words[3].split(":")[0]
		start = int(words[3].split(":")[1].split("-")[0])
		end = int(words[3].split(":")[1].split("-")[1])
		bin = get_bin_from_chr_loc(chr, start)

		hsh[bin] = [chr, start, end, bin]
		arr.append(hsh[bin])

	return hsh, arr

# given a trans-C ranking, return merfish loci closest to the top and bottom of the ranking
def get_top_bottom_clique_in_merfish(hsh_loci, clique):
	arr_in_merfish = []
	for b in clique:
		if b in hsh_loci:
			arr_in_merfish.append(b)
			
	return arr_in_merfish[0:40], arr_in_merfish[len(arr_in_merfish)-40:]

# for each trans-C ranking, find the proximity of the merfish loci at the top and bottom of the ranking,
# plot them and assess statistically the difference
def compute_top_v_bottom_proximity():
	im_med,im_fr,im_med_trans,im_med_cis,im_fr_trans,im_fr_cis,nlen = pickle_read(f"{T}mat_contact_IMR90_untreated.pkl")
	hsh_loci, arr_loci = get_loci_coordinates()
	weight_fr = get_proximity_weight(im_med, im_fr, top_loci, arr_loci)
	
	fo = open(f"{T}merfish_top_v_bottom.txt", "w")
	arr_top_weights, arr_bottom_weights = [], []

	for tf_id in get_tf_ids("IMR90"):
		m = read_tf_bed(tf_id)
		pi = read_pi(f"{P}Outs/TF_seed/{tf_id}.txt")
		clique = [x[0] for x in pi]
		
		# find the merfish loci at the top and bottom of the list
		top_loci, bottom_loci = get_top_bottom_clique_in_merfish(hsh_loci, clique)
		
		top_weight = get_proximity_weight(im_med, im_fr, top_loci, arr_loci)
		bottom_weight = get_proximity_weight(im_med, im_fr, bottom_loci, arr_loci)
		
		arr_top_weights.append(top_weight)
		arr_bottom_weights.append(bottom_weight)
		fo.write(f"{top_weight}\t{bottom_weight}\n")
	
	# are the two groups significantly different	
	p = stats.mannwhitneyu(arr_top_weights, tf_weight_med_arr, alternative="greater")
	print(f"p = {p}")


# compare the proximity of the trans-C cliques v random sets of loci in MERFISH data
def compare_transC_v_random_loci():
	im_med,im_fr,im_med_trans,im_med_cis,im_fr_trans,im_fr_cis,nlen = pickle_read(f"{T}mat_contact_IMR90_untreated.pkl")
	hsh_loci, arr_loci = get_loci_coordinates()
	weight_fr = get_proximity_weight(im_med, im_fr, top_loci, arr_loci)
	
	NUM_BINS = 40
	fo = open(f"{T}transc_v_random.txt", "w")
	
	arr_transC = []
	# for each trans-C clique
	for tf_id in get_tf_ids("IMR90"):
		m = read_tf_bed(tf_id)
		pi = read_pi(f"{P}Outs/TF_seed/{tf_id}_merfish.txt")
		clique = [x[0] for x in pi][0:NUM_BINS]
		
		weight = get_proximity_weight(im_med, im_fr, clique, arr_loci)
		arr_transC.append(weight)
		
		fo.write(f"{weight}\ttransC\n")

	
	# now draw 1000 times random loci and get their proximity
	arr_random = []
	for i in range(1000):
		rand_bins = list(np.random.choice(1041, NUM_BINS, replace=False))
		weight = get_proximity_weight(im_med, im_fr, rand_bins, arr_loci)
		arr_random.append(weight)
		
		fo.write(f"{weight}\trandom\n")


# find if cliques are significant in sprite data
def calculate_sprite_signifcant_cliques(a_hic, a_sprite):
	NUM_BINS = 40
	
	# find the matched random cliques in hic data
	TF_matched_random_seed(a_hic, 0.5, "GM12878", 1000)

	fo = open(f"{T}table_sprite.txt", "w")
	for tf_id in get_tf_ids("GM12878"):
		m = read_tf_bed(tf_id)
		pi = read_pi(f"{P}Outs/TF_seed/{tf_id}.txt")
		clique = [x[0] for x in pi][0:NUM_BINS]
		
		weight_clique = find_clique_weight(clique, a_sprite)
		
		# find all 1000 matched random weights and assess significance
		matched_random_weights = []
		for i in range(1000):
			pi = read_pi(f"{P}Outs/TF_seed/{tf_id}_mm{i}.txt")
			clique = [x[0] for x in pi][0:NUM_BINS]
		
			weight = find_clique_weight(clique, a_sprite)
			matched_random_weights.append(weight)
		
		pval = matched_random_significant(weight_clique, matched_random_weights)
		
		fo.write(f"{tf_id}\t{weight_clique}\t{np.avg(matched_random_weights)}\t{pval}\n")
		

# for each tf clique, assess if it is signinicant in SPRTIE data, Hi-C, both, or neither
def compare_hic_sprite_significant():
	hsh_sprite, hsh_hic = {}, {}
	for line in open(f"{T}table_hic.txt"):
		words = line.rstrip().split("\t")
		hsh_hic[words[0]] = [float(words[1]), float(words[2]), float(words[3])]
		
	for line in open(f"{T}table_sprite.txt"):
		words = line.rstrip().split("\t")
		hsh_sprite[words[0]] = [float(words[1]), float(words[2]), float(words[3])]
		
	fo = open(f"{T}sprite_hic_significant.txt", "w")
	for tf in hsh_hic:
		if hsh_hic[tf][2] < 0.05:
			if hsh_sprite[tf][2] < 0.05:
				fo.write(f"{hsh_hic[tf][0]/hsh_hic[tf][1]}\t{hsh_sprite[tf][0]/hsh_sprite[tf][1]}\tboth\n")
			else:
				fo.write(f"{hsh_hic[tf][0]/hsh_hic[tf][1]}\t{hsh_sprite[tf][0]/hsh_sprite[tf][1]}\thic_only\n")
		elif hsh_sprite[tf][2] < 0.05:
			fo.write(f"{hsh_hic[tf][0]/hsh_hic[tf][1]}\t{hsh_sprite[tf][0]/hsh_sprite[tf][1]}\tsprite_only\n")
		else:
			fo.write(f"{hsh_hic[tf][0]/hsh_hic[tf][1]}\t{hsh_sprite[tf][0]/hsh_sprite[tf][1]}\tneither\n")


# compare the weight of the trans-C cliques in SPRITE data to that of a random sets of loci
def compare_transC_v_random_loci_sprite(a_sprite):
	NUM_BINS = 40
	fo = open(f"{T}transc_v_random_sprite.txt", "w")
	
	arr_transC = []
	# for each trans-C clique
	for tf_id in get_tf_ids("GM12878"):
		pi = read_pi(f"{P}Outs/TF_seed/{tf_id}.txt")
		clique = [x[0] for x in pi][0:NUM_BINS]
		
		weight_hic = find_clique_weight(clique, a_hic)
		arr_transC.append(weight)
		
		fo.write(f"{weight}\ttransC\n")

	arr_random = []
	for i in range(1000):
		rand_bins = list(np.random.choice(nbins, NUM_BINS, replace=False))
		weight = find_clique_weight(rand_bins, a_sprite)
		arr_random.append(weight)
		
		fo.write(f"{weight}\trandom\n")


if __name__ == "__main__":
	pass
	
	
