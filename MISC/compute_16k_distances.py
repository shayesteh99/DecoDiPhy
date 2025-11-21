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
import re
from collections import deque
import gzip

def __label_tree__(tree_obj):
	is_labeled = True
	i = 0
	labels = set()
	for node in tree_obj.traverse_preorder():
		if node.is_leaf():
			continue
		if not node.label or node.label in labels or isinstance(node.label, float): 
			is_labeled = False
			node.label = 'I' + str(i)
			i += 1        
		labels.add(node.label)
	return is_labeled

def preprocess_tree(tree_obj, epsilon = 1e-4):
	# mean_edge = np.mean([n.edge_length for n in tree_obj.traverse_preorder() if not n.is_root()])
	# std = max(mean_edge/sequence_length, 1e-3)

	for n in tree_obj.traverse_preorder():
		if not n.is_root():
			if n.edge_length <= 0:
				n.edge_length = epsilon

def read_tree(tree):
	new_tree = re.sub(r"\{.*?\}", "", tree)
	tree_obj = read_tree_newick(new_tree)
	__label_tree__(tree_obj)
	preprocess_tree(tree_obj)
	return tree_obj

def merge_placements(placements, labels, is_weighted):
	branch_count = {}
	for p in placements:
		if len(p['p']) == 0:
			continue
		length = len(p['p'])
		total = length
		if is_weighted:
			total = 0
			for i in range(length):
				total += p['p'][i][4]
		for i in range(length):
			bid = p['p'][i][0]
			label = labels[bid]
			num = 1
			if is_weighted:
				num = p['p'][i][4]
			if label in branch_count:
				branch_count[label] += num / total
			else:
				branch_count[label] = num / total
	total = sum([branch_count[k] for k in branch_count])
	# print(total)
	branch_count = {k:branch_count[k]/total for k in branch_count}
	# print(branch_count)
	return branch_count

# def merge_placements(placements, is_weighted = False):
# 	branch_count = {}
# 	read_id = None
# 	for i in range(len(placements)):
# 		p = placements[i]
# 		# print(p[0], p[1])
# 		if p[1] == 'NaN':
# 			continue
# 		if not read_id:
# 			read_id = p[0]
# 			assgs = [p[1]]
# 		elif read_id == p[0]:
# 			assgs.append(p[1])
# 		else:
# 			# print(read_id, assgs)
# 			for a in assgs:
# 				if a in branch_count:
# 					branch_count[a] += 1 / len(assgs)
# 				else:
# 					branch_count[a] = 1 / len(assgs)
# 			read_id = p[0]
# 			assgs = [p[1]]

# 	for a in assgs:
# 		if a in branch_count:
# 			branch_count[a] += 1 / len(assgs)
# 		else:
# 			branch_count[a] = 1 / len(assgs)

# 	total = sum([branch_count[k] for k in branch_count])
# 	branch_count = {k:branch_count[k]/total for k in branch_count}

# 	print(sum([branch_count[k] for k in branch_count]))
# 	return branch_count

def get_terminal_branch_lengths(tree_obj):
	sum_lengths = {}
	leaf_counts = {}
	for n in tree_obj.traverse_postorder():
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

def place_queries(tree, branch_count):
	tree_obj = read_tree_newick(tree)
	terminal_lengths = get_terminal_branch_lengths(tree_obj) 
	# labels = [n.label for n in tree_obj.traverse_postorder()]
	# ids_to_labels = {i:labels[i] for i in range(len(labels))}

	# print(ids_to_labels)

	# edge_lengths = [n.edge_length for n in tree_obj.traverse_postorder() if not n.is_root()]
	# tlength = np.mean(edge_lengths)

	for q in branch_count:
		# label = ids_to_labels[q]
		node = tree_obj.label_to_node(selection='all')[q]
		tlength = terminal_lengths[q]
		if node.is_root():
			query = Node(label = 'q' + str(q), edge_length = tlength)
			parent = Node(label = 'p' + str(q), edge_length = 0)
			node.edge_length = 0
			parent.add_child(node)
			parent.add_child(query)
			tree_obj.root = parent
		else:
			query = Node(label = 'q' + str(q), edge_length = tlength)
			parent = Node(label = 'p' + str(q), edge_length = node.edge_length / 2)
			node.edge_length = node.edge_length / 2
			nparent = node.parent
			nparent.remove_child(node)
			nparent.add_child(parent)
			parent.add_child(node)
			parent.add_child(query)
	# print(tree_obj.newick())
	return tree_obj


def compute_avg_distances(tree_obj, branch_count):
	distance_matrix = tree_obj.distance_matrix(leaf_labels=True)
	avg_matrix = {}
	query_labels = ['q'+str(q) for q in branch_count]

	for q in branch_count:
		d = distance_matrix['q'+str(q)]
		for l in d:
			if l not in query_labels:
				if l not in avg_matrix:
					avg_matrix[l] = d[l] * branch_count[q]
				else:
					avg_matrix[l] += d[l] * branch_count[q]

	return avg_matrix


# def compute_avg_distances(tree_obj, branch_count):
# 	labels = [l.label for l in tree_obj.traverse_preorder()]

# 	distance_matrix = tree_obj.distance_matrix(leaf_labels=True)
# 	avg_matrix = {}
# 	query_labels = [q for q in branch_count if q in leaves]

# 	total_p = sum(branch_count[q] for q in query_labels)
# 	print(total_p)

# 	for q in query_labels:
# 		d = distance_matrix[q]
# 		for l in d:
# 			if l not in avg_matrix:
# 				avg_matrix[l] = d[l] * branch_count[q] / total_p
# 			else:
# 				avg_matrix[l] += d[l] * branch_count[q] / total_p

# 	return avg_matrix

def get_subtree(tree_obj, node, max_nodes = 500, max_bls = 2):
	q = deque()
	q.append((node, 0))

	num_nodes = 0
	sum_bls = 0
	visited = set()
	while num_nodes < max_nodes and sum_bls < max_bls and len(q) > 0:
		curr = q.popleft()
		num_nodes += 1
		sum_bls += curr[1]
		visited.add(curr[0])
		for c in curr[0].child_nodes():
			if c not in visited:
				q.append((c, c.edge_length))

		if curr[0].parent and curr[0].parent not in visited:
			q.append((curr[0].parent, curr[0].edge_length))

	# print(num_nodes)
	# print(sum_bls)
	# print([v.label for v in visited])

	for v in visited:
		for c in v.child_nodes():
			if c not in visited:
				v.remove_child(c)

		p = v.parent
		if p and p not in visited:
			p.remove_child(v)
			tree_obj.root = v

	# print(tree_obj.newick())
	

def match_trees(tree_obj, query_tree, branch_count):
	id_to_clade = {}

	num_to_id = {}

	ntree_leaves = [l.label for l in tree_obj.traverse_leaves()]
	qtree_leaves = [l.label for l in query_tree.traverse_leaves()]
	for n in tree_obj.traverse_preorder():
		n_leaves = [l.label for l in n.traverse_leaves()]
		# if len(n_leaves) < len(ntree_leaves) - len(n_leaves):
		id_to_clade[n.label] = set(n_leaves)
		# else:
		# 	id_to_clade[n.label] = set(ntree_leaves) - set(n_leaves)

	for n in query_tree.traverse_preorder():
		n_leaves = [l.label for l in n.traverse_leaves()]
		set_leaves = set(n_leaves)
		# if len(n_leaves) >= len(qtree_leaves) - len(n_leaves):
		# 	set_leaves = set(qtree_leaves) - set(n_leaves)

		is_l = False
		for k in id_to_clade:
			if id_to_clade[k] == set_leaves:
				num_to_id[k] = n.label
				is_l = True
				break 
		if not is_l:
			print(n.label, set_leaves)
		

	labels = [n.label for n in tree_obj.traverse_postorder()]
	ids_to_labels = {i:labels[i] for i in range(len(labels))}

	new_branch_count = {}

	root_id = None
	child_ids = []
	
	for i in range(len(labels)):
		if labels[i] == tree_obj.root.label:
			root_id = i
			continue

		for c in tree_obj.root.child_nodes():
			if labels[i] == c.label:
				child_ids.append(i)

	
	if root_id in branch_count:
		p = branch_count[root_id]
		for c in child_ids:
			if c not in branch_count:
				branch_count[c] = p/len(child_ids)
			else:
				branch_count[c] += p/len(child_ids)
		del branch_count[root_id]

	for b in branch_count:
		# print(b, ids_to_labels[b], num_to_id[ids_to_labels[b]])
		new_branch_count[num_to_id[ids_to_labels[b]]] = branch_count[b]

	return new_branch_count
	# print(new_branch_count)


def divide_tree(tree_obj, branch_count):
	sorted_branch_count = list(dict(sorted(branch_count.items(), key=lambda item: item[1], reverse=True)))

	n = len([n for n in tree_obj.traverse_preorder()])

	labels = [n.label for n in tree_obj.traverse_postorder()]
	ids_to_labels = {i:labels[i] for i in range(len(labels))}
	label_to_nodes = tree_obj.label_to_node(selection='all')

	visited_nodes = set()
	i = 0
	while len(visited_nodes) < n:
		node = label_to_nodes[ids_to_labels[sorted_branch_count[i]]]
		print(node.label)
		get_subtree(tree_obj, node, max_nodes = 10, max_bls = 2)
		i+=1
		return



def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input Placements")
	parser.add_argument('-t', '--tree', required=True, help="Input Tree")
	# parser.add_argument('-w', '--weighted', required=False, default='0', help="For multiplacements only")
	# parser.add_argument('-d', '--input_dir', required=True, help="Input dir")
	# parser.add_argument('-s', '--seed', required=True, help="Random Seed")
	# parser.add_argument('-s2', '--seed2', required=True, help="Random Seed 2")
	# parser.add_argument('-k', '--k', required=True, help="Number of query taxa")
	# parser.add_argument('-n', '--noise', required=True, help="Noise type")
	# parser.add_argument('-l', '--seq_length', required=False, default='1000', help="Sequence Length")
	# parser.add_argument('-w', '--weighted', required=False, default='0', help="For multiplacements only")
	# parser.add_argument('-c', '--correction', required=False, default='0', help="Correct the distances")
	# parser.add_argument('-f', '--fix_k', required=False, default='0', help="Fix k")
	# parser.add_argument('-m', '--method', required=False, choices=['hill', 'exhaustive', 'closest', 'closest_iterative'], default='hill', help="Seach method")
	# parser.add_argument('-p', '--placement_file', required=False, default='./.txt', help="placement file")
	parser.add_argument('-o', '--output', required=True, help="Output file")

	args = parser.parse_args()

	# with gzip.open(args.input, 'rt') as f:

	branch_count = {}
	with open(args.input,'r') as f:
		lines = f.readlines()
		lines = [l.split() for l in lines]

	for l in lines:
		branch_count[l[0]] = float(l[1])
	total = sum([branch_count[k] for k in branch_count])
	branch_count = {k:branch_count[k]/total for k in branch_count}

	with open(args.tree,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	qtree = read_tree_newick(inputTree)

	q_labels = [n.label for n in qtree.traverse_postorder()]
	q_branch_count = {b:branch_count[b] for b in branch_count if b in q_labels}
	total_p = sum([q_branch_count[b] for b in q_branch_count])

	print(total_p)
	# with open(os.path.join(outdir, "total_p.txt"),'w') as f:
	# 	f.write(str(total_p)+"\n")

	q_branch_count = {b:q_branch_count[b]/total_p for b in q_branch_count}

	tree_with_queries = place_queries(inputTree, q_branch_count)

	avg_matrix = compute_avg_distances(tree_with_queries, q_branch_count)

	with open(args.output,'w') as f:
		for l in avg_matrix:
			f.write(l + "\t" + str(avg_matrix[l]) + "\n")



	return 

	tree = placements['tree']
	tree_obj = read_tree(tree)

	is_weighted = args.weighted == '1'

	labels = [n.label for n in tree_obj.traverse_postorder()]
	branch_count = merge_placements(placements['placements'], labels, is_weighted)

	with open(args.output,'w') as f:
		for l in branch_count:
			f.write(l + "\t" + str(branch_count[l]) + "\n")

	return
	for t in range(16):
		outdir = os.path.join(args.outdir, "tree_" + str(t) + "_q1_r2")
		with open(os.path.join(outdir, "pruned_tree.trees"),'r') as f:
			inputTree = f.read().strip().split("\n")[0]

		qtree = read_tree_newick(inputTree)

		q_labels = [n.label for n in qtree.traverse_postorder()]
		q_branch_count = {b:branch_count[b] for b in branch_count if b in q_labels}
		total_p = sum([q_branch_count[b] for b in q_branch_count])

		print(total_p)
		with open(os.path.join(outdir, "total_p.txt"),'w') as f:
			f.write(str(total_p)+"\n")

		q_branch_count = {b:q_branch_count[b]/total_p for b in q_branch_count}

		tree_with_queries = place_queries(inputTree, q_branch_count)

		avg_matrix = compute_avg_distances(tree_with_queries, q_branch_count)

		with open(os.path.join(outdir, "distances.txt"),'w') as f:
			for l in avg_matrix:
				f.write(l + "\t" + str(avg_matrix[l]) + "\n")

	return
	
	avg_matrix = compute_avg_distances(tree_obj, branch_count)

	with open(args.outdir, 'w') as f:
		for l in avg_matrix:
			f.write(l + "\t" + str(avg_matrix[l]) + "\n")

	branch_count = {}
	with open(args.input,'r') as f:
		lines = f.readlines()
		lines = [l.split() for l in lines]
		for l in lines:
			branch_count[l[0]] = float(l[1])

	with open(args.input_tree,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	tree_obj = read_tree_newick(inputTree)

	avg_matrix = compute_avg_distances(tree_obj, branch_count)

	with open(args.outdir, 'w') as f:
		for l in avg_matrix:
			f.write(l + "\t" + str(avg_matrix[l]) + "\n")

	return



	tree = placements['tree']
	tree_obj = read_tree(tree)
	# print(tree_obj.newick())

	with open(args.input_tree,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	qtree = read_tree_newick(inputTree)

	is_weighted = args.weighted == '1'

	branch_count = merge_placements(placements['placements'], is_weighted)
	# print(branch_count)

	new_branch_count = match_trees(tree_obj, qtree, branch_count)

	return


	# print(dict(sorted(branch_count.items(), key=lambda item: item[1], reverse=True)))

	tree_with_queries = place_queries(tree_obj.newick(), branch_count)
	# print(tree_with_queries.newick())

	# divide_tree(tree_obj, branch_count)
	avg_matrix = compute_avg_distances(tree_with_queries, branch_count)

	with open(args.outdir, 'w') as f:
		for l in avg_matrix:
			f.write(l + "\t" + str(avg_matrix[l]) + "\n")


if __name__ == "__main__":
	main()  