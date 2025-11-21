import sys
import os
import argparse
import time
import numpy as np
from scipy.optimize import root
from scipy.optimize import least_squares
from scipy.optimize import minimize
from treeswift import *
import treeswift
import random
# import matplotlib.pyplot as plt
import cvxpy as cp
import json
from itertools import combinations
from collections import deque
from scipy.stats import chi2
from scipy.stats import wasserstein_distance
from scipy.spatial.distance import braycurtis 

def __label_tree__(tree_obj):
	is_labeled = True
	i = 0
	labels = set()
	for node in tree_obj.traverse_preorder():
		# if not node.is_root():
		# 	if node.edge_length < 0:
		# 		node.edge_length = 0
		if node.is_leaf():
			continue
		if not node.label or node.label in labels or node.label[0] != 'I': 
			is_labeled = False
			node.label = 'I' + str(i)
			i += 1        
		labels.add(node.label)
	return is_labeled

def preprocess_tree(tree_obj, sequence_length = 10**3):
	mean_edge = np.mean([n.edge_length for n in tree_obj.traverse_preorder() if not n.is_root()])
	std = max(mean_edge/sequence_length, 1e-3)

	for n in tree_obj.traverse_preorder():
		if not n.is_root():
			if n.edge_length < 1e-4:
				noise = np.random.normal(0, std)
				while noise < 1e-4:
					noise = np.random.normal(0, std)
				n.edge_length += noise
	if len(tree_obj.root.child_nodes()) == 2:
		sum_bls = sum(n.edge_length for n in tree_obj.root.child_nodes())
		for n in tree_obj.root.child_nodes():
			n.edge_length = sum_bls / 2

def distance_between(u, v):
        if u == v:
            return 0.
        elif u == v.parent:
            return v.edge_length
        elif v == u.parent:
            return u.edge_length
        u_dists = {u:0.}; v_dists = {v:0.}
        c = u; p = u.parent # u traversal
        while p is not None:
            u_dists[p] = u_dists[c]
            if c.edge_length is not None:
                u_dists[p] += c.edge_length
            if p == v:
            	return u_dists[p]
            c = p; p = p.parent
        c = v; p = v.parent # v traversal
        while p is not None:
            v_dists[p] = v_dists[c]
            if c.edge_length is not None:
                v_dists[p] += c.edge_length
            if p in u_dists:
                return u_dists[p] + v_dists[p]
            c = p; p = p.parent
        raise RuntimeError("u and v are not in the same Tree")

def get_random_query_taxa(tree_obj, leaves, k):
	taken = set()
	query = []
	q = None
	anchor = None

	while len(query) < k:
		while not q or q in taken or anchor in taken:
			q = random.choice(leaves)


		while q in taken or q.parent in taken or q.parent.parent in taken:
			q = random.choice(leaves)

		if q.parent.is_root():
			for c in q.parent.child_nodes():
				if c != q:
					sibling = c
					break
			if sibling.child_nodes()[0] in taken:
				continue
			else:
				query.append(q)
				taken.add(q)
				taken.add(sibling.child_nodes()[0])
		else:
			query.append(q)
			taken.add(q)
			taken.add(q.parent)
			if q.parent.parent:
				taken.add(q.parent.parent)

	return query

def get_neighbors(tree, node_to_index, radius = 2):
	node_neighbors = {}

	for n in tree.traverse_preorder():
		visited = set()
		q = deque()
		q.append((n, 0))
		visited.add(n.label)
		while len(q) > 0:
			curr = q.popleft()
			if curr[1] > radius:
				break

			for c in curr[0].child_nodes():
				if c.label not in visited:
					q.append((c, curr[1] + 1))
					visited.add(c.label)
			p = curr[0].parent
			if p and p.label not in visited:
				q.append((p, curr[1] + 1))
				visited.add(p.label)

			if p:
				for c in p.child_nodes():
					if c != curr[0] and c.label not in visited:
						q.append((c, curr[1] + 1))
						visited.add(c.label)

		visited.remove(n.label)
		node_neighbors[node_to_index[n.label]] = [node_to_index[v] for v in visited]
	return node_neighbors

def get_neighbors_for_node(tree,  node, radius):
	node_neighbors = {}

	visited = set()
	q = deque()
	q.append((node, 0))
	visited.add(node.label)
	node_neighbors[0] = [node.label]
	while len(q) > 0:
		curr = q.popleft()
		if curr[1] > radius:
			break

		for c in curr[0].child_nodes():
			if c.label not in visited:
				q.append((c, curr[1] + 1))
				if curr[1] + 1 not in node_neighbors:
					node_neighbors[curr[1] + 1] = []
				node_neighbors[curr[1] + 1].append(c.label)
				visited.add(c.label)
		p = curr[0].parent
		if p and p.label not in visited:
			q.append((p, curr[1] + 1))
			if curr[1] + 1 not in node_neighbors:
				node_neighbors[curr[1] + 1] = []
			node_neighbors[curr[1] + 1].append(p.label)
			visited.add(p.label)

		if p:
			for c in p.child_nodes():
				if c != curr[0] and c.label not in visited:
					q.append((c, curr[1] + 1))
					if curr[1] + 1 not in node_neighbors:
						node_neighbors[curr[1] + 1] = []
					node_neighbors[curr[1] + 1].append(c.label)
					visited.add(c.label)

	for n in node_neighbors:
		if tree.root.label in node_neighbors[n]:
			node_neighbors[n].remove(tree.root.label)
	return node_neighbors


def get_terminal_branch_lengths(tree_obj):
	sum_lengths = {}
	leaf_counts = {}
	for n in tree_obj.traverse_postorder():
		if n.is_root():
			n.edge_length = 0
		if n.is_leaf():
			sum_lengths[n.label] = 0
			leaf_counts[n.label] = 1
		else:
			sum_lengths[n.label] = 0
			leaf_counts[n.label] = 0
			for c in n.child_nodes():
				sum_lengths[n.label] += sum_lengths[c.label] + leaf_counts[c.label] * c.edge_length
				leaf_counts[n.label] += leaf_counts[c.label]

	terminal_lengths = {}
	for n in tree_obj.traverse_postorder():
		terminal_lengths[n.label] = sum_lengths[n.label] / leaf_counts[n.label] + n.edge_length / 2

	return terminal_lengths


def place_queries(tree, assignments):
	tree_obj = read_tree_newick(tree)
	terminal_lengths = get_terminal_branch_lengths(tree_obj) 
	labels = [n.label for n in tree_obj.traverse_postorder()]

	queries = {}
	labels_to_nodes = tree_obj.label_to_node(selection='all')
	for q in assignments:
		queries['q' + q] = assignments[q]
		node = labels_to_nodes[q]
		tlength = terminal_lengths[q]
		if node.is_root():
			query = Node(label = 'q' + q, edge_length = tlength)
			parent = Node(label = 'p' + q, edge_length = 0)
			node.edge_length = 0
			parent.add_child(node)
			parent.add_child(query)
			tree_obj.root = parent
		else:
			query = Node(label = 'q' + q, edge_length = tlength)
			parent = Node(label = 'p' + q, edge_length = node.edge_length / 2)
			node.edge_length = node.edge_length / 2
			nparent = node.parent
			nparent.remove_child(node)
			nparent.add_child(parent)
			parent.add_child(node)
			parent.add_child(query)

	return tree_obj, queries


def compute_avg_distances(tree_obj, queries):
	distance_matrix = tree_obj.distance_matrix(leaf_labels=True)
	avg_matrix = {}

	for q in queries:
		d = distance_matrix[q]
		for l in d:
			if l not in queries:
				if l not in avg_matrix:
					avg_matrix[l] = d[l] * queries[q]
				else:
					avg_matrix[l] += d[l] * queries[q]

	return avg_matrix


def compute_metrics(noisy_tree, noisy_queries, queries, input_tree, outdir):
	tree_obj = read_tree_newick(noisy_tree)
	for l in tree_obj.traverse_leaves():
		if l.label in noisy_queries or l.label in queries:
			l.edge_length = 0

	unifrac = compute_unifrac(tree_obj, queries, noisy_queries)
	weighted_unifrac = compute_weighted_unifrac(tree_obj, queries, noisy_queries)

	jaccard = jacard(set([q for q in queries]), set([q[1:] for q in noisy_queries]))

	true_abunds = {n.label:0 for n in input_tree.traverse_preorder()}
	noisy_abunds = {n.label:0 for n in input_tree.traverse_preorder()}

	for q in queries:
		true_abunds[q] = queries[q]
	true_abunds = [true_abunds[i] for i in true_abunds]

	for q in noisy_queries:
		noisy_abunds[q[1:]] = noisy_queries[q]
	noisy_abunds = [noisy_abunds[i] for i in noisy_abunds]

	wasser_dist = wasserstein_distance(true_abunds, noisy_abunds)
	bc_dist = braycurtis(true_abunds, noisy_abunds)

	estimated_k = len(noisy_queries)

	with open(os.path.join(outdir, "metrics.txt"), 'w') as f:
		f.write(str(jaccard) + "\t" + str(unifrac) + "\t" + str(weighted_unifrac) + "\t" + str(estimated_k) + "\t" + str(wasser_dist) + "\t" + str(bc_dist))

	# print(jaccard, "\t", unifrac, "\t", weighted_unifrac, "\t", estimated_k, "\t", wasser_dist, "\t",  bc_dist)

def get_diameter(tree):
	tree_obj = read_tree_newick(tree)
	for n in tree_obj.traverse_preorder():
		if not n.is_root():
			n.edge_length = 1
	return tree_obj.diameter()

def write_distances_to_file(distances, path):
	with open(path, 'w') as f:
		for a in distances:
			f.write(a +  "\t" + str(distances[a]) + "\n")

def write_queries_to_file(path, anchors, p, x = None, y = None):
	with open(path, 'w') as f:
		for i in range(len(anchors)):
			if x:
				f.write(anchors[i] + "\t" + str(p[i]) + "\t" + str(x[i]) + "\t" + str(y[i]) + "\n")
			else:
				f.write(anchors[i] + "\t" + str(p[i]) + "\n")

def save_tree_to_file(tree_newick, path):
	with open(path, 'w') as f:
		f.write(tree_newick)

def simulate_queries(tree_obj, init_queries, k, outdir, add_noise = 0, read_count = 10 ** 5, exp_scale = 1):
	leaves = [n for n in tree_obj.traverse_leaves()]

	if init_queries:
		query = [tree_obj.label_to_node(selection='all')[q[0]] for q in init_queries]
	else:
		query = get_random_query_taxa(tree_obj, leaves, k)
	query_labels = [q.label for q in query]

	index_to_leaf = [l.label for l in leaves if l not in query]
	leaf_to_index = {index_to_leaf[i]: i for i in range(len(index_to_leaf))}

	__label_tree__(tree_obj)
	# print(tree_obj.newick())
	# preprocess_tree(tree_obj)
	# print(tree_obj.newick())

	avg_matrix = {}

	#pick abundances
	if init_queries:
		true_p = [float(q[1]) for q in init_queries]
	else:
		true_p = np.random.rand(k)
		true_p /= np.sum(true_p)

		# while min(true_p) < 0.01/k:
		while min(true_p) < 0.1/(k**2):
			true_p = np.random.rand(k)
			true_p /= np.sum(true_p)

	queries = {query_labels[i]: true_p[i] for i in range(len(query_labels))}
	print(queries)

	with open(os.path.join(outdir, "query_labels.txt"), "w") as f:
		for q in queries:
			f.write(q + "\t" + str(queries[q]) + "\n")

	avg_matrix = compute_avg_distances(tree_obj, queries)
	# print(len(avg_matrix))
	# print(queries)

	#pruning the query taxa from the tree
	correct_anchor = []
	true_y = []
	true_x = []
	for q in query:
		parent = q.get_parent()
		parent.remove_child(q)
		sibling = parent.child_nodes()[0]
		if parent.is_root():
			true_y.append(q.edge_length + sibling.edge_length)
			true_x.append(1)
			parent.remove_child(sibling)
			tree_obj.root = sibling
			correct_anchor.append(sibling.child_nodes()[0].label)
		else:
			true_y.append(q.edge_length)
			true_x.append(sibling.edge_length / (sibling.edge_length + parent.edge_length))
			parent.contract()
			correct_anchor.append(sibling.label)

	# print(tree_obj.newick())
	write_queries_to_file(os.path.join(outdir, "true_queries.txt"), correct_anchor, true_p, true_x, true_y)
	save_tree_to_file(tree_obj.newick(), os.path.join(outdir, "pruned_tree.trees"))

	labels_to_nodes = tree_obj.label_to_node(selection='all')

	if add_noise == 1:
		assignments = {}
		for i in range(len(correct_anchor)):
			samples = np.random.exponential(scale=exp_scale * get_diameter(tree_obj.newick())/20 , size=int(read_count * true_p[i])).astype(int)
			# print(samples)
			if correct_anchor[i] not in assignments:
				assignments[correct_anchor[i]] = len(samples[samples == 0])
			else:
				assignments[correct_anchor[i]] += len(samples[samples == 0])
			radius = max(samples)
			node = labels_to_nodes[correct_anchor[i]]
			neighbors = get_neighbors_for_node(tree_obj, node, radius = radius)
			max_distance = max([a for a in neighbors])
			samples = samples[samples <= max_distance]
			for r in range(1, max(samples) + 1):
				count = len(samples[samples == r])
				reads = np.random.choice(neighbors[r], size=count, replace=True)
				for n in neighbors[r]:
					l = len(reads[reads == n])
					if l > 0:
						if n not in assignments:
							assignments[n] = l
						else:
							assignments[n] += l

		total_reads = sum([assignments[a] for a in assignments])
		print(total_reads)
		# for a in assignments:
		# 	print(a, "\t", assignments[a])
		assignments = {a:assignments[a]/total_reads for a in assignments}
		noisy_tree_obj, noisy_queries = place_queries(tree_obj.newick(), assignments)


		write_queries_to_file(os.path.join(outdir, "noisy_queries.txt"), [a for a in assignments], [assignments[a] for a in assignments])

		# compute_metrics(noisy_tree_obj.newick(), noisy_queries, queries, tree_obj, outdir)
		avg_matrix = compute_avg_distances(noisy_tree_obj, noisy_queries)

	write_distances_to_file(avg_matrix, os.path.join(outdir, "distances.txt"))


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input tree")
	parser.add_argument('-s', '--seed', required=True, help="Random Seed")
	parser.add_argument('-k', '--k', required=True, help="Number of query taxa")
	parser.add_argument('-q', '--queries', required=False, help="Query File")
	parser.add_argument('-n', '--noise', required=True, help="Noise type")
	parser.add_argument('-e', '--exp_scale', required=False, default='1', help="Exponential Scale")
	parser.add_argument('-o', '--outdir', required=False, default='', help="Output dir")

	args = parser.parse_args()

	random.seed(a=int(args.seed))
	np.random.seed(int(args.seed))

	start = time.time()

	with open(args.input,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	k = int(args.k)

	queries = None
	if args.queries:
		with open(args.queries,'r') as f:
			lines = f.readlines()
			queries = [l.split() for l in lines]

	simulate_queries(read_tree_newick(inputTree), queries, k, args.outdir, add_noise = int(args.noise), read_count = 10 ** 5, exp_scale = float(args.exp_scale))

	end_sim = time.time()
	print("Simulation Runtime: ", end_sim - start)

if __name__ == "__main__":
	main()  