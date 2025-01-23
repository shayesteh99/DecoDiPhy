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

def get_input_matrices(tree_obj, k, outdir, add_noise = 0, sequence_length = 10 ** 3):
	leaves = [n for n in tree_obj.traverse_leaves()]
	query = get_random_query_taxa(tree_obj, leaves, k)
	query_labels = [q.label for q in query]
	print(query_labels)
	index_to_leaf = [l.label for l in leaves if l not in query]
	leaf_to_index = {index_to_leaf[i]: i for i in range(len(index_to_leaf))}

	__label_tree__(tree_obj)
	# print(tree_obj.newick())
	# preprocess_tree(tree_obj)
	# print(tree_obj.newick())

	if add_noise == 1:
		noisy_tree_obj = read_tree_newick(tree_obj.newick())
		for n in noisy_tree_obj.traverse_preorder():
			if not n.is_root():
				noise = np.random.normal(0, n.edge_length/sequence_length)
				n.edge_length += noise

		dist_matrix = noisy_tree_obj.distance_matrix(leaf_labels=True)

	elif add_noise == 0:
		dist_matrix = tree_obj.distance_matrix(leaf_labels=True)

	avg_matrix = {}
	true_p = np.random.rand(k)
	true_p /= np.sum(true_p)

	while min(true_p) < 0.01/k:
		true_p = np.random.rand(k)
		true_p /= np.sum(true_p)


	for i in range(k):
		if add_noise == 2:
			noisy_tree_obj = read_tree_newick(tree_obj.newick())
			for n in noisy_tree_obj.traverse_preorder():
				if not n.is_root():
					n.edge_length += np.random.normal(0, n.edge_length/sequence_length)
			dist_matrix = noisy_tree_obj.distance_matrix(leaf_labels=True)
		q = query_labels[i]
		d = dist_matrix[q]
		for l in d:
			if l not in query_labels:
				if l not in avg_matrix:
					avg_matrix[l] = d[l] * true_p[i]
				else:
					avg_matrix[l] += d[l] * true_p[i]
	
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

	#matrix C
	index_to_node = [n for n in tree_obj.labels(leaves=True, internal=True)]
	node_to_index = {index_to_node[i]: i for i in range(len(index_to_node))}

	#save tree
	save_tree_to_file(tree_obj.newick(), node_to_index, outdir, "pruned_tree.trees")

	C = np.zeros((len(index_to_leaf), len(index_to_node)))

	for l in tree_obj.traverse_leaves():
		node = l
		while node:
			C[leaf_to_index[l.label], node_to_index[node.label]] = 1
			node = node.parent
	C[C == 0] = -1

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
	# print([k for k in new_dist_matrix])

	D = np.zeros((len(index_to_leaf), len(index_to_node)))
	for i in range(len(index_to_leaf)):
		for j in range(len(index_to_node)):
			leaf = index_to_leaf[i]
			node = index_to_node[j]
			D[i,j] = new_dist_matrix[node][leaf]

	#vector l
	l = np.array([label_to_node[i].edge_length for i in index_to_node])
	l[node_to_index[tree_obj.root.label]] = 1
	l = l[:,np.newaxis]

	#creating the input
	d = np.zeros(len(index_to_leaf))
	for i in range(len(index_to_leaf)):
		leaf = index_to_leaf[i]
		d[i] = avg_matrix[leaf]

	#get the true output
	# print(query_labels)
	# print(correct_anchor)
	# print(true_x)
	# print(true_y)
	# p = np.zeros(len(index_to_node))
	# x = np.zeros(len(index_to_node))
	# y = np.zeros(len(index_to_node))
	A = np.zeros((len(index_to_node), k))
	p = np.array(true_p)
	x = np.array(true_x)
	y = np.array(true_y)
	# print("true p: ", p)
	# print("true x: ", true_x)
	# print("true y: ", np.sum(true_y * true_p))

	# print(correct_anchor)
	# print(tree_obj.newick())
	for i in range(k):
		# print(query_labels[i])
		index = node_to_index[correct_anchor[i]]
		A[index,i] = 1
		# p[i] = 1/k

		# p[index] = 1/k
		# x[index] = true_x[i]
		# y[index] = true_y[i]
	

	# print(x)
	# np.savetxt("x.csv", x, delimiter=",")
	# term1 = D @ p
	# term2 = C @ (x * l * p)
	# term3 = np.ones(len(d)) * (np.sum(y) / k)
	# term3 = np.ones(len(d)) * np.sum(y * p)
	term1 = D @ A @ p
	# print(term1)
	term2 = C @ (A * l) @ (x * p)
	# print(term2)
	term3 = np.ones(len(d)) * np.sum(y * p)
	# print(term3)
	# Compute the residuals
	residuals = d - (term1 + term2 + term3)

	# print(np.sum(residuals**2))
	# print(residuals)
	# print(np.sum(residuals))

	# np.savetxt("D.csv", D, delimiter=",")
	# np.savetxt("input.csv", d, delimiter=",")
	# np.savetxt("l.csv", l, delimiter=",")
	# np.savetxt("C.csv", C, delimiter=",")


	return d, l, C, D, A, index_to_node, node_to_index, index_to_leaf, np.array([node_to_index[correct_anchor[i]] for i in range(k)]), np.array(true_x), np.sum(y * p), np.array(true_p), tree_obj

def solve_with_k_cvxpy(d, l, C, D, true_A, index_to_node, k, anchors, true_x, true_y, true_p):
	# print(A)
	# A[2][2] = 0
	# A[10][2] = 1
	n = len(index_to_node)
	L = len(d)

	# all_objectives = []
	# for _ in range(1000):
	# new_anchors = np.random.choice(n, k, replace = False)
	# if set(new_anchors) == set(anchors):
	# 	print(new_anchors, anchors)
	# 	continue
	# A = np.zeros((n,k))
	# for i in range(k):
	# 	A[new_anchors[i]][i] = 1
	# A = np.array(true_A)

	# p = cp.Variable(k)  
	p = np.array([1/k for _ in range(k)]) 
	# x = cp.Variable(k)
	x = np.array([1/(2*k) for _ in range(k)])         
	A = cp.Variable((n, k), boolean=True) 
	# y = cp.Variable()
	y = 0

	term1 = D @ A @ p   
	term2 = C @ (cp.multiply(A, l)) @ x
	term3 = np.ones(L) * y

	residuals = d - (term1 + term2 + term3)
	objective = cp.Minimize(cp.sum_squares(residuals))

	constraints = [
    # p >= 0,
    # x >= 0,
    # x <= p,
    # y >= 0,
  	# cp.sum(p) == 1,
  	cp.sum(A, axis=0) == 1,
  	cp.sum(A, axis=1) <= 1,
	]

	problem = cp.Problem(objective, constraints)
	problem.solve(solver=cp.GUROBI, verbose=True)

	# 	all_objectives.append(objective.value)
	print(objective.value)
	# print(min(all_objectives))

	# print(problem.status)

	# print("Optimal p:", p.value)
	# print("Optimal x:", x.value / p.value)
	print(anchors)
	print("Optimal A (binary):", np.where(A.value == 1)[0])
	# print(true_y)
	# print(y.value)

	# term1 = D @ np.array(A.value) @ p   
	# term2 = C @ (np.array(A.value) * l) @ (x * p)
	# term3 = np.ones(L) * y

	# residuals = d - (term1 + term2 + term3)
	# print(residuals)
	# print("Optimal y:", y.value)


def get_optimal_obj(d, l, C, D, anchors, index_to_node, k):
	t = time.time()
	DA = D[:, anchors]
	CAl = C[:, anchors] * l[anchors].T

	n = len(index_to_node)
	L = len(d)

	p = cp.Variable(k)   
	x = cp.Variable(k)   
	y = cp.Variable()

	term1 = DA @ p   
	term2 = CAl @ x
	term3 = np.ones(L) * y

	residuals = d - (term1 + term2 + term3)
	objective = cp.Minimize(cp.sum_squares(residuals))

	constraints = [
    p >= 0,
    x >= 0,
    x <= p,
    y >= 0,
  	cp.sum(p) == 1,
	]

	problem = cp.Problem(objective, constraints)
	try:
		problem.solve(verbose=False)
	except cp.error.SolverError:
		print("Solver reached max iterations. Using the last available solution.")
		return float('Inf'), None, None, None

	e = time.time()
	return objective.value, p.value, np.array([min(i, 1) for i in np.fabs(x.value/p.value)]), y.value, e-t


def hill_climbing(d, l, C, D, index_to_node, k, all_rounds, initial_anchors = None):
	opt_times = []
	n = len(index_to_node)
	L = len(d)
	if initial_anchors is None:
		anchors = np.random.choice([i for i in range(1, n)], k, replace = False)
	else:
		new = np.random.choice([i for i in range(1, n) if i not in initial_anchors], k-len(initial_anchors), replace = False)
		anchors = np.array(list(initial_anchors) + list(new))

	print(anchors)

	original, p, x, y, t = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)
	opt_times.append(t)
	min_p, min_x, min_y = p, x, y
	min_obj = original
	rounds = 0
	while True:
		init = time.time()
		rounds += 1
		print("=" * 200)
		print(rounds)
		min_obj = original
		# min_anchors = np.array(anchors)
		og_anchors = np.array(anchors)

		for ind in range(k):
			i = k-ind-1
			# new_anchors = np.array(anchors)
			min_val = 1
			max_val = n-1
			min_anchor = anchors[i]
			for v in range(min_val, max_val + 1):
				if v in anchors:
					continue
				# new_anchors[i] = v
				anchors[i] = v
				new_obj, new_p, new_x, new_y, t = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)
				opt_times.append(t)
				if min_obj > new_obj:
					min_obj = new_obj
					min_anchor = v
					# min_anchors[i] = v
					min_p = new_p
					min_x = new_x
					min_y = new_y
			anchors[i] = min_anchor
			# print(anchors)
		print(min_obj)
		end_round = time.time()
		round_info = {}
		round_info["k"] = k
		round_info["rounds"] = rounds
		round_info["loss"] = min_obj
		round_info["anchors"] = list(map(int, anchors))
		round_info["p"] = list(min_p)
		round_info["x"] = list(min_x)
		round_info["y"] = float(min_y)
		round_info["runtime"] = end_round - init
		round_info["opttime"] = np.mean(opt_times)
		all_rounds.append(round_info)
		if set(og_anchors) == set(anchors):
			return anchors, min_obj, min_p, min_x, min_y
		original = min_obj
		og_anchors = anchors

	return anchors, min_obj, min_p, min_x, min_y

def exhaustive_search(d, l, C, D, index_to_node, k):
	# anchors = np.array([39, 78, 88])
	# original, p, x, y = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)

	# return anchors, original, p, x, y

	n = len(index_to_node)
	all_anchors = list(combinations(range(1, n), k))

	anchors = np.random.choice([i for i in range(1, n)], k, replace = False)
	min_anchors = np.array(anchors)
	original, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)
	min_p, min_x, min_y = p, x, y
	min_obj = original

	for anchors in all_anchors:
		anchors = np.array([a for a in anchors])
		new_obj, new_p, new_x, new_y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)

		if min_obj > new_obj:
			min_obj = new_obj
			min_anchors = np.array(anchors)
			min_p = new_p
			min_x = new_x
			min_y = new_y

	return min_anchors, min_obj, min_p, min_x, min_y

def k_closest_leaves(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k):
	indices = np.argsort(d)[:k]
	leaves = np.array(index_to_leaf)[indices]
	anchors = np.array([node_to_index[l] for l in leaves])
	obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)

	return anchors, obj, p, x, y

def k_closest_leaves_iterative(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k):
	final_anchors = []
	for i in range(k):
		# num_a = k - i
		indices = np.argsort(d)[:k]
		leaves = np.array(index_to_leaf)[indices]
		anchors = np.array([node_to_index[l] for l in leaves])
		obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)
		index = np.argmax(p)
		while anchors[index] in final_anchors:
			p[np.argmax(p)] = -1
			index = np.argmax(p)

		DA = D[:, anchors[index]]
		CAl = C[:, anchors[index]] * l[anchors[index]]
		dist = (DA + CAl * x[index] + y)/k
		d = (d - dist) * k / (k-1)
		final_anchors.append(anchors[index])
	
	obj, p, x, y, _ = get_optimal_obj(d, l, C, D, final_anchors, index_to_node, k)
	return np.array(final_anchors), obj, p, x, y



def co_optimize(d, l, C, D, index_to_node, k):
	n = len(index_to_node)
	L = len(d)

	# anchors = np.random.choice(n, k, replace = False)
	# A = np.zeros((n,k))
	# for i in range(k):
	# 	A[anchors[i]][i] = 1
	p = np.array([1/k for _ in range(k)]) 
	x = np.array([1/(2*k) for _ in range(k)])  
	y = 0

	opt_obj = float("inf")
	rounds = 0
	while opt_obj > 1e-10:
		rounds += 1
		print("="*200)
		print(rounds)

		##optimaize A
		A = cp.Variable((n, k), boolean=True)

		term1 = D @ A @ p   
		term2 = C @ (cp.multiply(A, l)) @ x
		term3 = np.ones(L) * y

		residuals = d - (term1 + term2 + term3)
		objective = cp.Minimize(cp.sum_squares(residuals))

		constraints = [cp.sum(A, axis=0) == 1, 
						cp.sum(A, axis=1) <= 1]
		problem = cp.Problem(objective, constraints)
		problem.solve(solver = cp.GUROBI, verbose=False)

		print("second step:")
		print(objective.value)
		opt_obj = objective.value

		A = np.array(A.value)
		print(np.where(A == 1)[0])

		##optimize x,y,p
		p = cp.Variable(k)   
		x = cp.Variable(k)         
		y = cp.Variable()

		term1 = D @ A @ p   
		term2 = C @ (cp.multiply(A, l)) @ x
		term3 = np.ones(L) * y

		residuals = d - (term1 + term2 + term3)
		objective = cp.Minimize(cp.sum_squares(residuals))

		constraints = [
	    p >= 0,
	    x >= 0,
	    x <= p,
	    y >= 0,
	  	cp.sum(p) == 1,
		]

		problem = cp.Problem(objective, constraints)
		problem.solve(verbose=False)

		print("first step:")
		print(objective.value)

		y = np.array(y.value)
		x = np.array(x.value)
		p = np.array(p.value)

		print(p)
		print(x)
		print(y)

	return A, x, y, p


def solve_with_k(d, l, C, D, A, index_to_node, true_p, k):
	n = len(index_to_node)
	L = len(d)

	def objective(vars):
		# x = vars[:n]
		# y = vars[n:2*n]
		# p = vars[2*n:]  
		p = vars[:k]
		x = vars[k:2*k]
		# A = vars[2*k:2*k + (n*k)].reshape((n, k))
		# y = vars[2*k + (n*k)]
		y = vars[2*k]

		term1 = D @ A @ p
		term2 = C @ (A * l) @ (x * p)
		term3 = np.ones(L) * y

		residuals = d - (term1 + term2 + term3)
		return np.sum(residuals**2)

	def column_sum_constraints(vars):
	    A = vars[2*k:2*k + (n*k)].reshape((n, k))
	    col_sums = np.sum(A, axis=0)
	    return col_sums - 1

	# constraints = [
	# 				# {'type': 'eq', 'fun': lambda vars: np.sum(vars[:k]) - 1},
	# 				{'type': 'eq', 'fun': lambda vars: column_sum_constraints(vars)}]

	bounds = [(1/k, 1/k) for _ in range(k)] + [(0, 1) for _ in range(k,2*k)] + [(0, None)]
			# [(0, 1) for _ in range(2*k,2*k + (n*k))] +
	# bounds = [(0, 1) for _ in range(n)] + [(0, None) for _ in range(n, 2 * n)] + [(0, 0) for _ in range(2 * n, 3*n)]
	# for i in range(len(true_p)):
	# 	bounds[2*n + true_p[i]] = (1/3,1/3)
	# print(len(bounds))
	# initial_guesses = [np.abs(np.random.rand(3 * n)) for _ in range(10)]

	final_p = np.zeros(n)
	for _ in range(10):
		# g = np.abs(np.random.rand(2*k + (n*k) + 1))
		g = np.abs(np.random.rand(2*k + 1))
		print("+"*200)
		solution = minimize(objective, g, method='SLSQP', bounds=bounds)
			# , constraints=constraints)

		# Extract results
		p, x, y = solution.x[:k], solution.x[k:2*k], solution.x[2*k]#, solution.x[2*k + (n*k)]
		# A = solution.x[2*k:2*k + (n*k)].reshape((n, k))
		print("x:", x)
		print("y:", y)
		# print("P:", p)
		# print("A: ", A)
		print("Sum of P:", np.sum(p)) 
		print("Success:", solution.success)


		term1 = D @ A @ p
		term2 = C @ (A * l) @ (x * p)
		term3 = np.ones(L) * y

		# Compute the residuals
		residuals = d - (term1 + term2 + term3)
		print("Final residuals:", residuals)
		# print(d)
		print("Sum of squared residuals:", np.sum(residuals**2))

		# plt.plot(p, color="blue")
		# for i in range(len(true_p)):
		# 	plt.vlines(true_p[i], 0, 1/len(true_p), color="red")
		# plt.show()

		# final_p += p

	print("Final Results:")
	plt.plot(final_p / 10, color="blue")
	for i in range(len(true_p)):
		plt.vlines(true_p[i], 0, 1/len(true_p), color="red")
	plt.show()


def solve(d, l, C, D, index_to_node, true_p):
	n = len(index_to_node)
	L = len(d)

	def objective(vars):
		# x = vars[:n]
		# y = vars[n:2*n]
		# p = vars[2*n:]  
		p = vars[:n]
		x = vars[n:2*n]
		y = vars[2*n]

		term1 = D @ p
		term2 = C @ (x * l * p)
		term3 = np.ones(L) * y

		residuals = d - (term1 + term2 + term3)
		return np.sum(residuals**2)

	constraints = [{'type': 'eq', 'fun': lambda vars: np.sum(vars[:n]) - 1}]

	bounds = [(0, 1) for _ in range(n)] + [(0, 1) for _ in range(n)] + [(0, None)]
	# bounds = [(0, 1) for _ in range(n)] + [(0, None) for _ in range(n, 2 * n)] + [(0, 0) for _ in range(2 * n, 3*n)]
	# for i in range(len(true_p)):
	# 	bounds[2*n + true_p[i]] = (1/3,1/3)
	# print(len(bounds))
	# initial_guesses = [np.abs(np.random.rand(3 * n)) for _ in range(10)]

	final_p = np.zeros(n)
	for _ in range(10):
		g = np.abs(np.random.rand(2 * n + 1))
		print("+"*200)
		solution = minimize(objective, g, method='SLSQP', bounds=bounds, constraints=constraints)

		# Extract results
		# x, y, p = solution.x[:n], solution.x[n:2*n], solution.x[2*n:]
		p, x, y = solution.x[:n], solution.x[n:2*n], solution.x[2*n]
		print("x:", x[true_p])
		print("y:", y)
		print("P:", p)
		print("Sum of P:", np.sum(p)) 
		print("Success:", solution.success)


		term1 = D @ p
		term2 = C @ (x * l * p)
		term3 = np.ones(L) * y

		# Compute the residuals
		residuals = d - (term1 + term2 + term3)
		print("Final residuals:", residuals)
		# print(d)
		print("Sum of squared residuals:", np.sum(residuals**2))

		plt.plot(p, color="blue")
		for i in range(len(true_p)):
			plt.vlines(true_p[i], 0, 1/len(true_p), color="red")
		plt.show()

		final_p += p

	print("Final Results:")
	plt.plot(final_p / 10, color="blue")
	for i in range(len(true_p)):
		plt.vlines(true_p[i], 0, 1/len(true_p), color="red")
	plt.show()


def jacard(a, b):
	return len(a.intersection(b)) / len(a.union(b))

def find_maximum_curvature(all_obj):
	vals = np.log10(all_obj)
	if len(vals) < 3:
		if vals[0] - vals[-1] < 2:
			return 2
		return len(vals) + 1

	first = [vals[i] - vals[i+1] for i in range(len(vals) - 1)]
	second = [first[i] - first[i+1] for i in range(len(first) - 1)]
	max_curv = np.argmax(second)

	print(second)

	if np.max(second) < 0.1:
		print(vals[0] - vals[-1])
		if vals[0] - vals[-1] < 2:
			return 2
		return len(vals) + 1
	return max_curv + 3

def find_maximum_needle(all_obj):
	vals = np.log10(all_obj)
	m = (vals[-1] - vals[0]) / (len(vals) -1)
	b = vals[0]

	dist = [m*i - vals[i] + b for i in range(len(vals))]
	print(dist)

	length = np.sqrt((vals[-1] - vals[0])**2 + (len(vals) -1)**2)

	print(np.max(dist) / length)
	if np.max(dist) / length < 0.07:
		if vals[0] - vals[-1] <=1:
			return 2
		return len(vals) + 1

	max_needle = np.argmax(dist)
	return max_needle + 2


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input tree")
	parser.add_argument('-s', '--seed', required=True, help="Random Seed")
	# parser.add_argument('-s2', '--seed2', required=True, help="Random Seed 2")
	parser.add_argument('-k', '--k', required=True, help="Number of query taxa")
	parser.add_argument('-n', '--noise', required=True, help="Noise type")
	parser.add_argument('-l', '--seq_length', required=False, default='1000', help="Sequence Length")
	parser.add_argument('-f', '--fix_k', required=False, default='0', help="Fix k")
	parser.add_argument('-m', '--method', required=False, choices=['hill', 'exhaustive', 'closest', 'closest_iterative'], default='hill', help="Seach method")
	parser.add_argument('-o', '--outdir', required=False, default='1000', help="Output dir")

	args = parser.parse_args()

	random.seed(a=int(args.seed))
	np.random.seed(int(args.seed))

	start = time.time()

	with open(args.input,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	k = int(args.k)

	sequence_length	= 10**3
	if args.seq_length:
		sequence_length = int(args.seq_length)
	d, l, C, D, A, index_to_node, node_to_index, index_to_leaf, anchors, true_x, true_y, true_p, pruned_tree = get_input_matrices(read_tree_newick(inputTree), k, args.outdir, add_noise = int(args.noise), sequence_length = sequence_length)
	sorted_indices = np.argsort(anchors)
	anchors = anchors[sorted_indices]
	true_p = true_p[sorted_indices]
	true_x = true_x[sorted_indices]
	print("anchors: ", anchors)
	print("true_p: ", true_p)
	print("true_x: ", true_x)
	print("true_y: ", true_y)

	end_sim = time.time()
	print("Simulation Runtime: ", end_sim - start)

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

	if args.fix_k == '1':
		if args.method == "exhaustive":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = exhaustive_search(d, l, C, D, index_to_node, k)
		elif args.method == "closest":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = k_closest_leaves(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k)
		elif args.method == "closest_iterative":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = k_closest_leaves_iterative(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k)
		elif args.method == "hill":
			opt_anchors = None
			for i in range(2, k+1):
				opt_anchors, opt_obj, opt_p, opt_x, opt_y = hill_climbing(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors)
		print("opt_anchors: ", opt_anchors)
		print("opt_x: ", opt_x)
		print("opt_p: ", opt_p)
		print("opt_y: ", opt_y)

	elif args.fix_k == '0':
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
			if min(opt_p) < 1/len(index_to_node) and i > 2:
				print("stopped by p0")
				break

			all_obj.append(opt_obj)
			all_x.append(opt_x)
			all_p.append(opt_p)
			all_y.append(opt_y)
			all_anchors.append(opt_anchors)
			if opt_obj < 1e-10:
				print("stopped by obj0")
				break

			if prev_obj and (prev_obj - opt_obj)/prev_obj < 0.1:
				print("stopped by obj_imp")
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
				print("stopped by obj1/2")
				break

	# if np.fabs(max(opt_x) - 1) < 1e-3:
	# 	label = index_to_node[opt_anchors[np.argmax(opt_x)]]
	# 	node = pruned_tree.label_to_node(selection='all')[label]
	# 	parent = node_to_index[node.parent.label]
	# 	if not node.parent.is_root():
	# 		sisters = []
	# 		for c in node.parent.child_nodes():
	# 			if c != node:
	# 				sisters.append(node_to_index[c.label])
	# 		sisters.append(parent)
	# 		for s in sisters:
	# 			anchor = opt_anchors[np.argmax(opt_x)]
	# 			if s in opt_anchors:
	# 				if args.fix_k == '1':
	# 					continue
	# 				new_anchors = np.array([a for a in opt_anchors if a != anchor])
	# 			else:
	# 				new_anchors = np.array(opt_anchors)
	# 				new_anchors[np.argmax(opt_x)] = s
	# 			new_obj, new_p, new_x, new_y = get_optimal_obj(d, l, C, D, new_anchors, index_to_node, len(new_anchors))
	# 			if np.fabs(max(new_x) - 1) > 1e-3:
	# 				opt_anchors = new_anchors
	# 				opt_x = new_x
	# 				opt_y = new_y
	# 				opt_p = new_p
	# 				break

	if max(opt_x) > 1 - 1e-3:
		indices = np.where(opt_x > 1 - 1e-3)[0]
		for i in indices:
			label = index_to_node[opt_anchors[i]]
			node = pruned_tree.label_to_node(selection='all')[label]
			parent = node_to_index[node.parent.label]
			if not node.parent.is_root():
				sisters = []
				for c in node.parent.child_nodes():
					if c != node:
						sisters.append(node_to_index[c.label])
				sisters.append(parent)
				for s in sisters:
					if s in opt_anchors:
						continue
					new_anchors = np.array(opt_anchors)
					new_anchors[i] = s
					new_obj, new_p, new_x, new_y, _ = get_optimal_obj(d, l, C, D, new_anchors, index_to_node, len(new_anchors))
					if new_x[i] < 1 - 1e-3:
						opt_anchors = new_anchors
						opt_x = new_x
						opt_y = new_y
						opt_p = new_p
						break

	if min(opt_x) < 1e-4:
		indices = np.where(opt_x < 1e-4)[0]
		for i in indices:
			label = index_to_node[opt_anchors[i]]
			node = pruned_tree.label_to_node(selection='all')[label]
			for c in node.child_nodes():
				s = node_to_index[c.label]
				if s in opt_anchors:
					continue
				new_anchors = np.array(opt_anchors)
				new_anchors[i] = s
				new_obj, new_p, new_x, new_y, _ = get_optimal_obj(d, l, C, D, new_anchors, index_to_node, len(new_anchors))
				if new_x[i] > 1e-4 and np.fabs(new_x[i] - 1) > 1e-3:
					opt_anchors = new_anchors
					opt_x = new_x
					opt_y = new_y
					opt_p = new_p
					break

	root_edges = [node_to_index[c.label] for c in pruned_tree.root.child_nodes()]
	if len(root_edges) == 2 and root_edges[0] in opt_anchors and root_edges[1] in opt_anchors:
		if args.fix_k == '0':
			list(opt_anchors).remove(root_edges[1])
			opt_anchors = np.array(opt_anchors)
			opt_obj, opt_p, opt_x, opt_y, _ = get_optimal_obj(d, l, C, D, opt_anchors, index_to_node, len(opt_anchors))

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