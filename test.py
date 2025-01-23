from treeswift import *
import numpy as np

def preprocess_tree(tree_obj, epsilon = 1e-4):
	for n in tree_obj.traverse_preorder():
		if not n.is_root():
			if n.edge_length <= 0:
				n.edge_length = epsilon

with open("../scaffolds_backbone.newick",'r') as f:
	inputTree = f.read().strip().split("\n")[0]


count = 13
tree = read_tree_newick(inputTree)
subtree_roots = []
for n in tree.traverse_preorder():
	if n.is_root() or n.is_leaf():
		continue

	is_sub = False
	ancestors = [a for a in n.traverse_ancestors()]
	for a in ancestors:
		if a in subtree_roots:
			is_sub = True
			break
	if is_sub:
		continue

	subtree = tree.extract_subtree(n)
	leaves = [l for l in subtree.labels(internal = False)]
	if len(leaves) >= 30:
		d = subtree.diameter()
		if d < 0.3:
			edge_lengths = np.array([n.edge_length for n in subtree.traverse_preorder() if not n.is_root()])
			ratio = len(edge_lengths[edge_lengths < 0.001]) / len(edge_lengths)
			if ratio < 0.5:
				print(len(leaves), d)
				print(subtree.newick())
				print("=" * 200)
				# print(ancestors)
				# print(subtree_roots)
				count += 1
				subtree_roots.append(n)
				preprocess_tree(subtree)
				# subtree.write_tree_newick("./10k_dataset/tree_" + str(count) + ".tre")


exit()
dist_matrix = {}

with open("test.txt",'r') as f:
	lines = f.readlines()
	# inputTree = f.read().strip().split("\n")[0]
	lines = [l.split() for l in lines]
	# for l in lines[0]:
	# 	dist_matrix[l] = {}
	for line in lines[1:]:
		genome = line[0]
		for i in range(len(line[1:])):
			if float(line[i+1]) < 0.25 and genome != lines[0][i]:
				if genome not in dist_matrix:
					dist_matrix[genome] = {}
				dist_matrix[genome][lines[0][i]] = float(line[i+1])

for d in dist_matrix:
	print(d , dist_matrix[d])

# for l in dist_matrix:
# 	for g in dist_matrix[l]:

