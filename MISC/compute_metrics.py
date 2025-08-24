import sys
import os
import argparse
import time
from treeswift import *
import treeswift
import json
import numpy as np
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

def place_queries(tree_obj, original_tree, anchors, p, x, prefix):
	query_labels = {}

	for i in range(len(anchors)):
		anchor = tree_obj.label_to_node(selection='all')[str(anchors[i])]
		origianl_anchor = original_tree.label_to_node(selection='all')[str(anchors[i])]
		l = origianl_anchor.edge_length

		leaf = Node(label = prefix + str(i), edge_length = 0)
		query_labels[prefix + str(i)] = p[i]

		if anchor.parent.label == origianl_anchor.parent.label:
			parent = anchor.parent
			new_parent = Node(label = "P" + prefix + str(i), edge_length = (1 - x[i]) * l)
			parent.remove_child(anchor)
			parent.add_child(new_parent)

			new_parent.add_child(leaf)
			new_parent.add_child(anchor)
			anchor.edge_length = x[i] * l
			continue

		if anchor.edge_length > x[i] * l:
			parent = anchor.parent
			new_parent = Node(label = "P" + prefix + str(i), edge_length = anchor.edge_length - x[i] * l)
			parent.remove_child(anchor)
			parent.add_child(new_parent)

			new_parent.add_child(leaf)
			new_parent.add_child(anchor)
			anchor.edge_length = x[i] * l

		else:
			new_anchor = anchor.parent
			length = new_anchor.edge_length
			parent = new_anchor.parent
			new_parent = Node(label = "P" + prefix + str(i), edge_length = 0)
			parent.remove_child(new_anchor)
			parent.add_child(new_parent)

			new_parent.add_child(leaf)
			new_parent.add_child(new_anchor)
			new_anchor.edge_length = x[i] * l - anchor.edge_length
			new_parent.edge_length = length - new_anchor.edge_length

	return query_labels

def prune_non_query_leaves(tree_obj, true_labels, final_labels):
	for l in tree_obj.traverse_leaves():
		if not l.label in true_labels and not l.label in final_labels:
			parent = l.parent
			if parent.is_root():
				parent.remove_child(l)
				sibling = parent.child_nodes()[0]
				parent.remove_child(sibling)
				tree_obj.root = sibling
			else:
				parent.remove_child(l)
				parent.contract()
	print(tree_obj.newick())

def compute_weighted_unifrac(tree_obj, true_labels, final_labels):
	true_abunds = {}
	final_abunds = {}

	u = 0
	D = 0
	for n in tree_obj.traverse_postorder():
		if n.is_root():
			break
		if n.is_leaf():
			true_abunds[n.label] = 0
			final_abunds[n.label] = 0
			if n.label in true_labels:
				true_abunds[n.label] = true_labels[n.label]
			if n.label in final_labels:
				final_abunds[n.label] = final_labels[n.label]	
		else:
			true_abunds[n.label] = 0
			final_abunds[n.label] = 0
			for c in n.child_nodes():
				true_abunds[n.label] += true_abunds[c.label]
				final_abunds[n.label] += final_abunds[c.label]
		u += n.edge_length * np.fabs(true_abunds[n.label] - final_abunds[n.label])
		D += n.edge_length * (true_abunds[n.label] + final_abunds[n.label])
		# D += n.edge_length
	# print(u)
	return u/D

def compute_unifrac(tree_obj, true_labels, final_labels):
	true_abunds = {}
	final_abunds = {}

	u = 0
	D = 0
	for n in tree_obj.traverse_postorder():
		if n.is_root():
			break
		if n.is_leaf():
			true_abunds[n.label] = False
			final_abunds[n.label] = False
			if n.label in true_labels:
				true_abunds[n.label] = True
			if n.label in final_labels:
				final_abunds[n.label] = True
				
		else:
			true_abunds[n.label] = False
			final_abunds[n.label] = False
			for c in n.child_nodes():
				if true_abunds[c.label]:
					true_abunds[n.label] = True
				if final_abunds[c.label]:
					final_abunds[n.label] = True

		if true_abunds[n.label] and not final_abunds[n.label]:
			u += n.edge_length
		if not true_abunds[n.label] and final_abunds[n.label]:
			u += n.edge_length
		D += n.edge_length
	return u/D

def jacard(a, b):
	return len(a.intersection(b)) / len(a.union(b))

def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--tree', required=True, help="Input tree")
	parser.add_argument('-q', '--queries', required=True, help="True queries")
	parser.add_argument('-e', '--est_queries', required=True, help="Estimated queries")
	# parser.add_argument('-f', '--format', required=False, choices=['json', 'jplace', 'txt'], default='json', help="Estimated queries format")

	args = parser.parse_args()

	with open(args.tree,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	tree_obj = read_tree_newick(inputTree)

	root_edges = []
	if len(tree_obj.root.child_nodes()) == 2:
		root_edges = [c.label for c in tree_obj.root.child_nodes()]

	true_x = None
	true_anchors, true_p = [], []
	with open(args.queries, 'r') as f:
		lines = f.readlines()
		lines = [l.split() for l in lines]
		if len(lines[0]) > 2:
			true_x = []
		for l in lines:
			true_anchors.append(l[0])
			true_p.append(float(l[1]))
			if len(l) > 2:
				true_x.append(float(l[2]))
		if not true_x:
			true_x = [0.5 for _ in range(len(true_anchors))]

	true_labels = place_queries(tree_obj, read_tree_newick(inputTree), true_anchors, true_p, true_x, prefix = "T")

	est_format = args.est_queries.split(".")[-1]

	if est_format == "json":
		with open(args.est_queries, "r") as f:
			all_rounds = json.load(f)
			final_round = all_rounds[-1]
			est_anchors = final_round["anchors"]
			est_p = final_round["p"]
			est_x = final_round["x"]

	elif est_format == "txt":
		est_x = None
		est_anchors, est_p = [], []
		with open(args.est_queries, "r") as f:
			lines = f.readlines()
			lines = [l.split() for l in lines]
			if len(lines[0]) > 2:
				est_x = []
			for l in lines:
				est_anchors.append(l[0])
				est_p.append(float(l[1]))
				if len(l) > 2:
					est_x.append(float(l[2]))
		if not est_x:
			est_x = [0.5 for _ in range(len(est_anchors))]

	final_labels = place_queries(tree_obj, read_tree_newick(inputTree), est_anchors, est_p, est_x, prefix = "F")

	weighted_unifrac = compute_weighted_unifrac(tree_obj, true_labels, final_labels)
	unifrac = compute_unifrac(tree_obj, true_labels, final_labels)

	estimated_k = len(est_anchors)

	if len(root_edges) > 0:
		for i in range(len(true_anchors)):
			if true_anchors[i] == root_edges[1]:
				true_anchors[i] = root_edges[0]
		for i in range(len(est_anchors)):
			if est_anchors[i] == root_edges[1]:
				est_anchors[i] = root_edges[0]

	true_anchors = np.array(true_anchors)
	final_anchors = np.array(est_anchors)

	jaccard = jacard(set(true_anchors), set(final_anchors))

	true_abunds = {n.label:0 for n in read_tree_newick(inputTree).traverse_preorder()}
	final_abunds = {n.label:0 for n in read_tree_newick(inputTree).traverse_preorder()}

	for i in range(len(true_anchors)):
		true_abunds[true_anchors[i]] = true_p[i]
	true_abunds = [true_abunds[i] for i in true_abunds]

	for i in range(len(final_anchors)):
		final_abunds[str(final_anchors[i])] = est_p[i]
	final_abunds = [final_abunds[i] for i in final_abunds]


	wasser_dist = wasserstein_distance(true_abunds, final_abunds)


	bc_dist = braycurtis(true_abunds, final_abunds)
	print(jaccard, "\t", unifrac, "\t", weighted_unifrac, "\t", estimated_k, "\t", wasser_dist, "\t",  bc_dist)



if __name__ == "__main__":
	main()