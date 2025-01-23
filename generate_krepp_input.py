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

def save_tree_to_file(tree_newick, node_to_index, outdir, filename):
	tree_obj = read_tree_newick(tree_newick)
	for n in tree_obj.traverse_preorder():
		n.label = node_to_index[n.label]
	with open(os.path.join(outdir, filename), 'w') as f:
		f.write(tree_obj.newick())

def get_input_matrices(tree_obj, k, outdir, p_dist = "uniform"):
	leaves = [n for n in tree_obj.traverse_leaves()]
	query = get_random_query_taxa(tree_obj, leaves, k)
	query_labels = [q.label for q in query]
	# print(query_labels)
	index_to_leaf = [l.label for l in leaves if l not in query]
	leaf_to_index = {index_to_leaf[i]: i for i in range(len(index_to_leaf))}

	__label_tree__(tree_obj)
	# print(tree_obj.newick())
	# preprocess_tree(tree_obj)
	# print(tree_obj.newick())

	dist_matrix = tree_obj.distance_matrix(leaf_labels=True)

	avg_matrix = {}
	if p_dist == "uniform":
		true_p = np.random.rand(k)
		true_p /= np.sum(true_p)

		while min(true_p) < 0.01/k:
			true_p = np.random.rand(k)
			true_p /= np.sum(true_p)

	elif p_dist == "exp":
		true_p = np.random.exponential(size = k)
		true_p /= np.sum(true_p)

		while min(true_p) < 0.01/k:
			true_p = np.random.exponential(size = k)
			true_p /= np.sum(true_p)


	for i in range(k):
		q = query_labels[i]
		d = dist_matrix[q]
		for l in d:
			if l not in query_labels:
				if l not in avg_matrix:
					avg_matrix[l] = d[l] * true_p[i]
				else:
					avg_matrix[l] += d[l] * true_p[i]
	
	#pruning the query taxa from the tree
	# correct_anchor = []
	# true_y = []
	# true_x = []
	# for q in query:
	# 	parent = q.get_parent()
	# 	parent.remove_child(q)
	# 	sibling = parent.child_nodes()[0]
	# 	if parent.is_root():
	# 		true_y.append(q.edge_length + sibling.edge_length)
	# 		true_x.append(1)
	# 		parent.remove_child(sibling)
	# 		tree_obj.root = sibling
	# 		correct_anchor.append(sibling.child_nodes()[0].label)
	# 	else:
	# 		true_y.append(q.edge_length)
	# 		true_x.append(sibling.edge_length / (sibling.edge_length + parent.edge_length))
	# 		parent.contract()
	# 		correct_anchor.append(sibling.label)
	# print(tree_obj.newick())

	#matrix C
	index_to_node = [n for n in tree_obj.labels(leaves=True, internal=True)]
	node_to_index = {index_to_node[i]: i for i in range(len(index_to_node))}

	#save tree
	save_tree_to_file(tree_obj.newick(), node_to_index, outdir, "pruned_tree.trees")

	# C = np.zeros((len(index_to_leaf), len(index_to_node)))

	# for l in tree_obj.traverse_leaves():
	# 	node = l
	# 	while node:
	# 		C[leaf_to_index[l.label], node_to_index[node.label]] = 1
	# 		node = node.parent
	# C[C == 0] = -1

	#matrix D	
	new_dist_matrix = tree_obj.distance_matrix(leaf_labels=True)
	for l in tree_obj.traverse_leaves():
		dist = 0
		node = l
		while node:
			if node.label not in new_dist_matrix:
				new_dist_matrix[node.label] = {}
			new_dist_matrix[node.label][l.label] = dist
			if node.edge_length:
				dist += node.edge_length
			node = node.parent

	for i in range(len(index_to_leaf)):
		for j in range(len(index_to_node)):
			leaf = index_to_leaf[i]
			node = index_to_node[j]
			if leaf not in new_dist_matrix[node]:
				for l in new_dist_matrix[node]:
					if l in index_to_leaf:
						new_dist_matrix[node][leaf] = new_dist_matrix[leaf][l] - new_dist_matrix[node][l]
						break

	label_to_node = tree_obj.label_to_node(selection='all')


	#creating the input
	d = {}
	for i in range(len(index_to_leaf)):
		leaf = index_to_leaf[i]
		d[leaf] = avg_matrix[leaf]


	return d, query_labels, index_to_node, node_to_index, np.array(true_p), tree_obj


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input tree")
	parser.add_argument('-s', '--seed', required=True, help="Random Seed")
	# parser.add_argument('-s2', '--seed2', required=True, help="Random Seed 2")
	parser.add_argument('-k', '--k', required=True, help="Number of query taxa")
	# parser.add_argument('-n', '--noise', required=True, help="Noise type")
	# parser.add_argument('-l', '--seq_length', required=False, default='1000', help="Sequence Length")
	parser.add_argument('-o', '--outdir', required=False, default='1000', help="Output dir")
	parser.add_argument('-p', '--pdist', required=False, default='uniform', help="p Distribution")

	args = parser.parse_args()

	random.seed(a=int(args.seed))
	np.random.seed(int(args.seed))

	start = time.time()

	with open(args.input,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	k = int(args.k)

	# sequence_length	= 10**3
	# if args.seq_length:
	# 	sequence_length = int(args.seq_length)
	d, query_labels, index_to_node, node_to_index, true_p, pruned_tree = get_input_matrices(read_tree_newick(inputTree), k, args.outdir, p_dist = args.pdist)
	# sorted_indices = np.argsort(anchors)
	# anchors = anchors[sorted_indices]
	# true_p = true_p[sorted_indices]
	# true_x = true_x[sorted_indices]
	# print("anchors: ", anchors)
	# print("queries: ", query_labels)
	# print("true_p: ", true_p)
	# print(d)

	with open(os.path.join(args.outdir, "queries.txt"), 'w') as f:
		for i in range(len(query_labels)):
			f.write(query_labels[i] + "\t" + str(true_p[i]) + "\n")

	with open(os.path.join(args.outdir, "distances.txt"), 'w') as f:
		for l in d:
			f.write(l + "\t" + str(d[l]) + "\n")
	# print("true_x: ", true_x)
	# print("true_y: ", true_y)

	end_sim = time.time()
	print("Simulation Runtime: ", end_sim - start)
	# random.seed(a=int(args.seed2))
	# np.random.seed(int(args.seed2))

	return

	all_rounds = []
	round_info = {}
	round_info["k"] = k
	round_info["rounds"] = "true"
	round_info["loss"] = 0
	round_info["anchors"] = list(map(int, anchors))
	round_info["p"] = list(true_p)
	round_info["x"] = list(true_x)
	round_info["y"] = float(true_y)
	round_info["runtime"] = 0
	all_rounds.append(round_info)


	opt_anchors = None
	prev_obj = None
	all_obj = []
	all_x = []
	all_p = []
	all_y = []
	all_anchors = []
	for i in range(2, len(index_to_node)):
		print("+" * 200)
		print(i)
		opt_anchors, opt_obj, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors)
		print("opt_x: ", opt_x)
		print("opt_p: ", opt_p)
		print("opt_y: ", opt_y)
		if min(opt_p) < 1e-2 and i > 2:
			break

		all_obj.append(opt_obj)
		all_x.append(opt_x)
		all_p.append(opt_p)
		all_y.append(opt_y)
		all_anchors.append(opt_anchors)
		if opt_obj < 1e-10:
			break
		prev_obj = opt_obj

	print(all_obj)

	opt_anchors = all_anchors[-1]
	opt_x = all_x[-1]
	opt_y = all_y[-1]
	opt_p = all_p[-1]
	for i in range(len(all_obj)-1):
		if all_obj[i] / all_obj[-1] < 2:
			opt_anchors = all_anchors[i]
			opt_x = all_x[i]
			opt_y = all_y[i]
			opt_p = all_p[i]
			break

	if max(opt_x) > 1 or np.fabs(max(opt_x) - 1) < 1e-3:
		label = index_to_node[opt_anchors[np.argmax(opt_x)]]
		node = pruned_tree.label_to_node(selection='all')[label]
		parent = node_to_index[node.parent.label]
		if not node.parent.is_root():
			sisters = []
			for c in node.parent.child_nodes():
				if c != node:
					sisters.append(node_to_index[c.label])
			sisters.append(parent)
			for s in sisters:
				anchor = opt_anchors[np.argmax(opt_x)]
				if s in opt_anchors:
					opt_anchors = np.array([a for a in opt_anchors if a != anchor])
				else:	
					opt_anchors[np.argmax(opt_x)] = s
				opt_obj, opt_p, opt_x, opt_y = get_optimal_obj(d, l, C, D, opt_anchors, index_to_node, len(opt_anchors))
				if np.fabs(max(opt_x) - 1) > 1e-3:
					break

	end = time.time()

	sorted_indices = np.argsort(opt_anchors)
	opt_anchors = opt_anchors[sorted_indices]
	opt_p = opt_p[sorted_indices]
	opt_x = opt_x[sorted_indices]

	round_info = {}
	round_info["k"] = len(opt_anchors)
	round_info["rounds"] = "final"
	round_info["loss"] = all_rounds[-1]["loss"]
	round_info["anchors"] = list(map(int, opt_anchors))
	round_info["p"] = list(opt_p)
	round_info["x"] = list(opt_x)
	round_info["y"] = float(opt_y)
	round_info["runtime"] = end - end_sim
	all_rounds.append(round_info)

	print("optimal anchors: ", opt_anchors)
	print("optimal p: ", opt_p)
	print("optimal x: ", opt_x)
	print("optimal y: ", opt_y)

	print("Jaccard: ", jacard(set(anchors), set(opt_anchors)))

	print("Optimization Runtime: ", end - end_sim) 

	with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
		json.dump(all_rounds, f)

	return
	max_n = find_maximum_needle(all_obj)
	# max_n += 2
	print("true k:", k)
	print("ratio k: ", len(opt_anchors))
	# print("imp k: ", imp_k)
	# print("max curve k: ", max_c)
	print("max needle k: ", max_n)
	#draw all abj
	plt.plot([i+2 for i in range(len(all_obj))], [np.sum(np.log(p))/len(p) for p in all_p], color = "lime")
	plt.plot([i+2 for i in range(len(all_obj))], np.log10(all_obj))
	plt.vlines(k, min(np.log10(all_obj)), max(np.log10(all_obj)), color="grey")
	# plt.vlines(len(opt_anchors), min(np.log10(all_obj)), max(np.log10(all_obj)), color="purple")
	# plt.vlines(imp_k, min(np.log10(all_obj)), max(np.log10(all_obj)), color="pink")
	# plt.vlines(max_c, min(np.log10(all_obj)), max(np.log10(all_obj)), color="blue")
	# plt.vlines(max_n, min(np.log10(all_obj)), max(np.log10(all_obj)), color="cyan")
	# plt.plot([2, len(all_obj)+1], [np.log10(all_obj[0]), np.log10(all_obj[-1])], color = 'red')

	plt.show()
	# plt.savefig("fig.png")
	return

	# opt_anchors = None
	# for i in range(2, len(index_to_node)):
	# 	print("+" * 200)
	# 	print(i)
	# 	opt_anchors, res, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors)
	# 	print("opt_x: ", opt_x)
	# 	print("opt_p: ", opt_p)
	# 	if res:
	# 		print("loss is zero")
	# 		if np.fabs(max(opt_x) - 1) < 1e-6:
	# 			print("x is one.")
	# 			print(opt_x)
	# 			label = index_to_node[opt_anchors[np.argmax(opt_x)]]
	# 			node = pruned_tree.label_to_node(selection='all')[label]
	# 			parent = node_to_index[node.parent.label]
	# 			if node.parent.is_root():
	# 				break
	# 			sisters = []
	# 			for c in node.parent.child_nodes():
	# 				if c != node:
	# 					sisters.append(node_to_index[c.label])
	# 			sisters.append(parent)
	# 			for s in sisters:
	# 				opt_anchors[np.argmax(opt_x)] = s
	# 				print("+" * 200)
	# 				print(i)
	# 				opt_anchors, res, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors)
	# 				if np.fabs(max(opt_x) - 1) > 1e-6:
	# 					break
	# 		break
	# 	if min(opt_p) < 1e-3 and i > 2:
	# 		print("p is small")
	# 		print(opt_p)
	# 		print("decrease k")
	# 		index = np.argmin(opt_p)
	# 		opt_anchors, res, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i-1, all_rounds, initial_anchors = np.delete(opt_anchors, index))
	# 		break
		# if min(opt_x) < 1e-6:
		# 	print("x is small")
		# 	print(opt_x)
		# 	print("try children")
		# 	index = np.where(opt_x < 1e-6)
		# 	for ind in index[0]:
		# 		label = index_to_node[ind]
		# 		node = pruned_tree.label_to_node(selection='all')[label]
		# 		for c in node.child_nodes():
		# 			opt_anchors[ind] = node_to_index[c.label]
		# 			opt_anchors, res, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors)
		# 			if min(opt_x) > 1e-6:
		# 				print("increase k")
		# 				continue

		# if min(opt_x) < 1e-6 and i > 2:
		# 	index = np.argmin(opt_x)
		# 	if opt_p[index] < 1e-2:
		# 		print("x and p are small")
		# 		print(opt_p)
		# 		print(opt_x)
		# 		print("decrease k")
		# 		opt_anchors, res, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i-1, all_rounds, initial_anchors = np.delete(opt_anchors, index))
		# 		break
		# print("increase k")

	end = time.time()
	sorted_indices = np.argsort(opt_anchors)
	opt_anchors = opt_anchors[sorted_indices]
	opt_p = opt_p[sorted_indices]
	opt_x = opt_x[sorted_indices]

	round_info = {}
	round_info["k"] = len(opt_anchors)
	round_info["rounds"] = "final"
	round_info["loss"] = all_rounds[-1]["loss"]
	round_info["anchors"] = list(map(int, opt_anchors))
	round_info["p"] = list(opt_p)
	round_info["x"] = list(opt_x)
	round_info["y"] = float(opt_y)
	round_info["runtime"] = end - end_sim
	all_rounds.append(round_info)

	print("optimal anchors: ", opt_anchors)
	print("optimal p: ", opt_p)
	print("optimal x: ", opt_x)
	print("optimal y: ", opt_y)

	with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
		json.dump(all_rounds, f)

	print("Jaccard: ", jacard(set(anchors), set(opt_anchors)))

	print("Optimization Runtime: ", end - end_sim) 







if __name__ == "__main__":
	main()  