import sys
import os
import argparse
import time
from treeswift import *
import treeswift
import json
import numpy as np

def compute_p_and_r(true, est):
	precision = 0
	for i in est:
		if i in true:
			precision += 1
	recall = 0
	for i in true:
		if i in est:
			recall += 1
	return precision/len(est), recall/len(true)

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
	# print(u)
	return u/D

def place_queries(tree_obj, original_tree, qinfo, prefix):
	k = int(qinfo["k"])
	anchors = qinfo["anchors"]
	x = qinfo["x"]
	p = qinfo["p"]

	query_labels = {}

	for i in range(k):
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

def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input dir")

	args = parser.parse_args()

	#compute precision and recall

	with open(os.path.join(args.input, "pruned_tree.trees"),'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	# tree_obj = read_tree_newick(inputTree)

	with open(os.path.join(args.input, "all_rounds.json"), "r") as f:
		all_rounds = json.load(f)

	true_anchors = all_rounds[0]["anchors"]
	true_round = all_rounds[0]

	all_pr = {}

	for i in range(1, len(all_rounds)):
		k = all_rounds[i]["k"]
		if i == len(all_rounds) - 1:
			p, r = compute_p_and_r(true_anchors, all_rounds[i]["anchors"])
			tree_obj = read_tree_newick(inputTree)
			true_labels = place_queries(tree_obj, read_tree_newick(inputTree), true_round, prefix = "T")
			final_labels = place_queries(tree_obj, read_tree_newick(inputTree), all_rounds[i], prefix = "F")
			weighted_unifrac = compute_weighted_unifrac(tree_obj, true_labels, final_labels)
			pro = np.array(all_rounds[i]['p']) 
			all_pr[all_rounds[i]["k"]] = (p, r, weighted_unifrac, all_rounds[i]['loss'], len(pro[pro<0.001]))
			break

		next_k = all_rounds[i+1]["k"]
		if k != next_k:
			p, r = compute_p_and_r(true_anchors, all_rounds[i]["anchors"])
			tree_obj = read_tree_newick(inputTree)
			true_labels = place_queries(tree_obj, read_tree_newick(inputTree), true_round, prefix = "T")
			final_labels = place_queries(tree_obj, read_tree_newick(inputTree), all_rounds[i], prefix = "F")
			weighted_unifrac = compute_weighted_unifrac(tree_obj, true_labels, final_labels)
			pro = np.array(all_rounds[i]['p']) 
			all_pr[all_rounds[i]["k"]] = (p, r, weighted_unifrac, all_rounds[i]['loss'], len(pro[pro<0.001]))
	# print(all_pr)
	for k in all_pr:
		print(k, "\t", all_pr[k][0], "\t", all_pr[k][1], "\t", all_pr[k][2], "\t", all_pr[k][3], "\t", all_pr[k][4])
	return

	#compute running time

	with open(os.path.join(args.input, "pruned_tree.trees"),'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	tree_obj = read_tree_newick(inputTree)
	n_leaves = len([l for l in tree_obj.traverse_leaves()])

	with open(os.path.join(args.input, "all_rounds.json"), "r") as f:
		all_rounds = json.load(f)

	final_round = all_rounds[-1]
	total_runtime = final_round['runtime']
	estimated_k = len(final_round['anchors'])

	one_round = all_rounds[-2]
	one_round_runtime = one_round['runtime']
	opt_runtime = one_round['opttime']

	num_rounds = len(all_rounds) -2

	print(n_leaves, "\t", estimated_k, "\t", total_runtime, "\t", one_round_runtime, "\t", opt_runtime, "\t", num_rounds)

if __name__ == "__main__":
	main()