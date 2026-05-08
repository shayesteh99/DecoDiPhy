import sys
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import time
import numpy as np
from treeswift import *
import treeswift
import random
import cvxpy as cp
import osqp
import scipy.sparse as sp
import json
from itertools import combinations
from collections import deque
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from .jutil import extended_newick
from .__init__ import __version__
from .optimize import AnchorOSQP
import tracemalloc

print("DecoDiPhy", __version__)


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
		if not node.label or node.label in labels or isinstance(node.label, float): 
			is_labeled = False
			node.label = 'I' + str(i)
			i += 1        
		labels.add(node.label)
	return is_labeled

def preprocess_input(tree_obj, distances):
	total = len([l for l in tree_obj.traverse_leaves()])
	leaves = [l.label for l in tree_obj.traverse_leaves() if l.label in distances]
	for l in tree_obj.traverse_leaves():
		if l.label not in leaves:
			parent = l.parent
			parent.remove_child(l)
			if parent.is_root():
				child = parent.child_nodes()[0]
				parent.remove_child(child)
				child.edge_length = 0
				tree_obj.root = child
			else:
				parent.contract()
	
	distances = {d:distances[d] for d in distances if d in leaves}

	if len(leaves) < total:
		return True
	return False

def post_process_output(tree_obj, og_tree, opt_anchors, opt_x):
	new_anchors = []
	new_x = []
	for i in range(len(opt_anchors)):
		anchor = opt_anchors[i]
		x = opt_x[i]
		node = tree_obj.label_to_node(selection = "all")[anchor]
		og_node = og_tree.label_to_node(selection = "all")[anchor]
		if node.parent.label == og_node.parent.label:
			new_anchors.append(anchor)
			new_x.append(x)
		else:
			length = node.edge_length * x
			while length > og_node.edge_length:
				length -= og_node.edge_length
				og_node = og_node.parent
			new_anchors.append(og_node.label)
			new_x.append(length / og_node.edge_length)

	return new_anchors, new_x

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


def get_neighbors(tree, node_to_index, radius = 2):
	node_neighbors = {}

	for n in tree.traverse_preorder():
		if n.is_root():
			continue
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
		if tree.root.label in visited:
			visited.remove(tree.root.label)
		node_neighbors[node_to_index[n.label]] = [node_to_index[v] for v in visited]
	return node_neighbors

def get_immediate_neighbors(tree, node_to_index):
	node_neighbors = {}

	for n in tree.traverse_preorder():
		if n.is_root():
			continue
		node_neighbors[node_to_index[n.label]] = []
		for c in n.child_nodes():
			node_neighbors[node_to_index[n.label]].append(node_to_index[c.label])
		parent = n.parent
		if not parent.is_root():
			node_neighbors[node_to_index[n.label]].append(node_to_index[parent.label])
			for c in parent.child_nodes():
				if c != n:
					node_neighbors[node_to_index[n.label]].append(node_to_index[c.label])
		else:
			for c in parent.child_nodes():
				if c != n:
					for gc in c.child_nodes():
						node_neighbors[node_to_index[n.label]].append(node_to_index[gc.label])
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

def get_input_matrices(tree_obj, distances):
	leaves = [n for n in tree_obj.traverse_leaves()]
	index_to_leaf = [l.label for l in leaves]
	leaf_to_index = {index_to_leaf[i]: i for i in range(len(index_to_leaf))}

	index_to_node = [n for n in tree_obj.labels(leaves=True, internal=True) if n != tree_obj.root.label]
	node_to_index = {index_to_node[i]: i for i in range(len(index_to_node))}

	#matrix C
	C = np.zeros((len(index_to_leaf), len(index_to_node)))

	for l in tree_obj.traverse_leaves():
		node = l
		while not node.is_root():
			C[leaf_to_index[l.label], node_to_index[node.label]] = 1
			node = node.parent
	C[C == 0] = -1

	#matrix D	
	dist_matrix = tree_obj.distance_matrix(leaf_labels=True)
	for l in tree_obj.traverse_leaves():
		dist = 0
		node = l
		while node:
			if node.label not in dist_matrix:
				dist_matrix[node.label] = {}
			dist_matrix[node.label][l.label] = dist
			if node.edge_length:
				dist += node.edge_length
			node = node.parent

	for i in range(len(index_to_leaf)):
		for j in range(len(index_to_node)):
			leaf = index_to_leaf[i]
			node = index_to_node[j]
			if leaf not in dist_matrix[node]:
				for l in dist_matrix[node]:
					if l in index_to_leaf:
						dist_matrix[node][leaf] = dist_matrix[leaf][l] - dist_matrix[node][l]
						break

	label_to_node = tree_obj.label_to_node(selection='all')

	D = np.zeros((len(index_to_leaf), len(index_to_node)))
	for i in range(len(index_to_leaf)):
		for j in range(len(index_to_node)):
			leaf = index_to_leaf[i]
			node = index_to_node[j]
			D[i,j] = dist_matrix[node][leaf]

	#vector l
	l = np.array([label_to_node[i].edge_length for i in index_to_node])
	# l[node_to_index[tree_obj.root.label]] = 1
	l = l[:,np.newaxis]

	#creating the input
	d = np.zeros(len(index_to_leaf))
	for i in range(len(index_to_leaf)):
		leaf = index_to_leaf[i]
		d[i] = distances[leaf]

	return d, l, C, D, index_to_node, node_to_index, index_to_leaf

def build_cvxpy_problem(L, k):
    p = cp.Variable(k)
    x = cp.Variable(k)
    y = cp.Variable()

    DA_param = cp.Parameter((L, k))
    CAl_param = cp.Parameter((L, k))
    d_param = cp.Parameter(L)

    residuals = d_param - (DA_param @ p + CAl_param @ x + y)

    objective = cp.Minimize(cp.sum_squares(residuals))

    constraints = [
        p >= 0,
        x >= 0,
        x <= p,
        y >= 0,
        cp.sum(p) == 1,
    ]

    problem = cp.Problem(objective, constraints)

    res= {}
    res['problem'] = problem
    res['p'] = p
    res['x'] = x
    res['y'] = y
    res['DA'] = DA_param
    res['CAl'] = CAl_param
    res['d'] = d_param

    return res


def get_optimal_obj(problem, d, l_T, C_T, D_T, anchors, index_to_node, k):
	t = time.time()
	# DA = D[:, anchors]
	DA = D_T[anchors].T
	# CAl = C[:, anchors] * l[anchors].T
	CAl = C_T[anchors].T * l_T[anchors]

	n = len(index_to_node)
	L = len(d)

	# p = cp.Variable(k)   
	# x = cp.Variable(k)   
	# y = cp.Variable()

	problem['DA'].value = DA
	problem['CAl'].value = CAl
	problem['d'].value = d

	# term1 = DA @ p
	# term2 = CAl @ x
	# term3 = np.ones(L) * y

	# residuals = d - (term1 + term2 + term3)
	# objective = cp.Minimize(cp.sum_squares(residuals))

	# constraints = [
    # p >= 0,
    # x >= 0,
    # x <= p,
    # y >= 0,
  	# cp.sum(p) == 1,
	# ]

	# problem = cp.Problem(objective, constraints)
	# print(problem.size_metrics)
	# print(problem.size_metrics.__dict__)
	try:
		problem['problem'].solve(verbose=False)
	except cp.error.SolverError:
		print("Solver reached max iterations. Using the last available solution.")
		e = time.time()
		return float('Inf'), None, None, None, e-t

	e = time.time()
	return problem['problem'].value, problem['p'].value, np.array([min(i, 1) for i in np.fabs(problem['x'].value/problem['p'].value)]), problem['y'].value, e-t


def get_optimal_obj_osqp(d, l_T, C_T, D_T, anchors, index_to_node, k):
	t = time.time()
	# DA = D[:, anchors]
	DA = D_T[anchors].T
	# CAl = C[:, anchors] * l[anchors].T
	CAl = C_T[anchors].T * l_T[anchors]

	n = len(index_to_node)
	L = len(d)

	M = np.hstack([DA, CAl, np.ones((L,1))])
	P = 2.0 * (M.T @ M)
	# P = P + 1e-9 * np.eye(P.shape[0])
	q = -2.0 * (M.T @ d) 

	n_vars = 2*k + 1
	row_data = []
	col_data = []
	val_data = []
	l_list = []
	u_list = []
	r = 0

    # --- p >= 0 ---
	for i in range(k):
		row_data.append(r); col_data.append(i); val_data.append(1.0)
		l_list.append(0.0); u_list.append(np.inf)
		r += 1

	# --- x >= 0 ---
	for i in range(k):
		row_data.append(r); col_data.append(k + i); val_data.append(1.0)
		l_list.append(0.0); u_list.append(np.inf)
		r += 1

	# --- x <= p  -> x - p <= 0 ---
	for i in range(k):
		# -p + x <= 0  (i.e., lower=-inf, upper=0)
		row_data.append(r); col_data.append(i); val_data.append(-1.0)
		row_data.append(r); col_data.append(k + i); val_data.append(1.0)
		l_list.append(-np.inf); u_list.append(0.0)
		r += 1

	# --- y >= 0 ---
	row_data.append(r); col_data.append(2*k); val_data.append(1.0)
	l_list.append(0.0); u_list.append(np.inf)
	r += 1

	# --- sum(p) = 1 ---
	for i in range(k):
		row_data.append(r); col_data.append(i); val_data.append(1.0)
	l_list.append(1.0); u_list.append(1.0)
	r += 1

	A = sp.csc_matrix((val_data, (row_data, col_data)), shape=(r, n_vars))
	l_vec = np.array(l_list, dtype=float)
	u_vec = np.array(u_list, dtype=float)

	P_csc = sp.csc_matrix(P)

	prob = osqp.OSQP()
	prob.setup(P=P_csc, q=q, A=A, l=l_vec, u=u_vec, eps_abs=1e-6, eps_rel=1e-6, polish=True, verbose=False)
	res = prob.solve()
    
	if res.info.status_val not in [1, 2]:
		print("OSQP did not converge")
		return float('inf'), None, None, None, time.time()-t

	z_opt = res.x
	p_opt = z_opt[:k]
	x_opt = z_opt[k:2*k]
	y_opt = z_opt[-1]

	e = time.time()

	obj_val = np.sum((d - (DA @ p_opt + CAl @ x_opt + y_opt))**2)

	return obj_val, p_opt, np.array([min(i, 1) for i in np.fabs(x_opt/p_opt)]), y_opt, e - t



def hill_climbing(d, l_T, C_T, D_T, index_to_node, k, all_rounds, solver, anchors = None, quick='0', neighbors = None, p_thresh = 0.01):

	def build_DA_CAl(anchors):
		DA = D_T[anchors].T
		CAl = C_T[anchors].T * l_T[anchors]
		return DA, CAl


	opt_times = []
	opt_memos = []
	n = len(index_to_node)
	L = len(d)

	print(anchors)
	original, p, x, y, t, m = solver.solve(anchors)
	DA, CAl = build_DA_CAl(anchors)

	opt_times.append(t)
	opt_memos.append(m)
	min_p, min_x, min_y = p, x, y
	min_obj = original
	last_anchor = anchors[-1]
	rounds = 0
	while True:
		init = time.time()
		rounds += 1
		print("=" * 200)
		print(rounds)
		min_obj = original

		og_anchors = np.array(anchors)

		for ind in range(k):
			i = k-ind-1
			min_val = 0
			max_val = n-1
			min_anchor = anchors[i]

			##possible anchor choices
			choices = range(min_val, max_val + 1)
			if ind > 0 and k > 2 and quick == '1':
				choices = neighbors[anchors[i]]
			count = 0
			for v in choices:
				if v in anchors:
					continue

				anchors[i] = v


				new_obj, new_p, new_x, new_y, t, _ = solver.solve(anchors)
				count += 1

				opt_times.append(t)
				if min_obj > new_obj:
					min_obj = new_obj
					min_anchor = v
					min_p = new_p
					min_x = new_x
					min_y = new_y
					last_anchor = v


			anchors[i] = min_anchor
			new_obj, new_p, new_x, new_y, t, _ = solver.solve(anchors)

		print(min_obj, np.log10(min_obj))
		end_round = time.time()
		round_info = {}
		round_info["k"] = k
		round_info["rounds"] = rounds
		round_info["loss"] = min_obj
		round_info["anchors"] = [index_to_node[i] for i in anchors]
		round_info["p"] = list(min_p)
		round_info["x"] = list(min_x)
		round_info["y"] = float(min_y)
		round_info["runtime"] = end_round - init
		round_info["opttime"] = np.mean(opt_times)
		round_info["memory"] = np.mean(opt_memos)
		all_rounds.append(round_info)
		if set(og_anchors) == set(anchors):
			return anchors, min_obj, min_p, min_x, min_y, round_info
		original = min_obj
		og_anchors = anchors

	return anchors, min_obj, min_p, min_x, min_y, round_info

def reverse_hill_climbing(d, l, C, D, index_to_node):
	n = len(index_to_node)
	anchors = np.array([i for i in range(n)])
	k = len(anchors)
	obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)

	# while np.min(p) < 1/n:
	# 	anchors = anchors[p >= 1/n]
	# 	k = len(anchors)
	# 	obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)

	return anchors, obj, p, x, y

def reverse_hill_climbing_fixed_k(d, l, C, D, index_to_node, k, min_p=0.01):
	n = len(index_to_node)
	anchors = np.array([i for i in range(n)])
	obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, len(anchors))
	while np.min(p) < min_p:
		if len(anchors[p >= min_p]) < k:
			break
		anchors = anchors[p >= min_p]
		obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, len(anchors))

	while len(anchors) > k:
		i = np.argmin(p)
		anchors = np.delete(anchors, i)
		obj, p, x, y, _ = get_optimal_obj(d, l, C, D, anchors, index_to_node, len(anchors))

	return anchors, obj, p, x, y


def exhaustive_search(d, l, C, D, index_to_node, k):
	# anchors = np.array([39, 78, 88])
	# original, p, x, y = get_optimal_obj(d, l, C, D, anchors, index_to_node, k)

	# return anchors, original, p, x, y

	n = len(index_to_node)
	all_anchors = list(combinations(range(n), k))

	anchors = np.random.choice([i for i in range(n)], k, replace = False)
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

def jacard(a, b):
	return len(a.intersection(b)) / len(a.union(b))


def get_dist_matrix_new(tree_obj, node_to_index):
	new_tree = read_tree_newick(tree_obj.newick())
	for n in new_tree.traverse_preorder():
		if n.is_root():
			continue
		n.label = node_to_index[n.label]
		if n.parent.is_root():
			n.edge_length = 0.5
		else:
			n.edge_length = 1

	dist_matrix = new_tree.distance_matrix(leaf_labels = True)

	for l in new_tree.traverse_preorder():
		if l.is_root():
			continue
		if l.label not in dist_matrix:
			dist_matrix[l.label] = {}
		dist = 1
		node = l.parent
		while node and not node.is_root():
			if node.label not in dist_matrix:
				dist_matrix[node.label] = {}
			dist_matrix[node.label][l.label] = dist
			dist_matrix[l.label][node.label] = dist
			if node.edge_length:
				dist += node.edge_length
			node = node.parent

	for n in new_tree.traverse_internal():
		if n.is_root():
			continue
		for m in new_tree.traverse_internal():
			if n == m:
				# dist_matrix[m.label][n.label] = 0
				continue
			if m.is_root():
				continue

			if n.label in dist_matrix[m.label]:
				continue
			leaf_n = [l for l in n.traverse_leaves()][0]
			leaf_m = [l for l in m.traverse_leaves()][0]
			dist_matrix[m.label][n.label] = dist_matrix[leaf_m.label][leaf_n.label] - dist_matrix[leaf_m.label][m.label] - dist_matrix[leaf_n.label][n.label]
			dist_matrix[n.label][m.label] = dist_matrix[leaf_m.label][leaf_n.label] - dist_matrix[leaf_m.label][m.label] - dist_matrix[leaf_n.label][n.label]


	for n in new_tree.traverse_internal():
		if n.is_root():
			continue
		for m in new_tree.traverse_leaves():

			if n.label in dist_matrix[m.label]:
				continue

			leaf_n = [l for l in n.traverse_leaves()][0]
			dist_matrix[m.label][n.label] = dist_matrix[m.label][leaf_n.label] - dist_matrix[leaf_n.label][n.label]
			dist_matrix[n.label][m.label] = dist_matrix[m.label][leaf_n.label] - dist_matrix[leaf_n.label][n.label]

	return dist_matrix

def get_dist_matrix_patristic(tree_obj):
	dist_matrix = {}
	node_set = {}
	height = {}
	LCAs = {}
	for n in tree_obj.traverse_postorder():
		if n.is_leaf():
			node_set[n.label] = [n]
			height[n.label] = n.edge_length
			continue

		if n.is_root():
			node_set[n.label] = [n]
			for c1 in range(len(n.children) -1):
				c1_nodes = node_set[n.children[c1].label]
				for c2 in range(c1+1,len(n.children)):
					c2_nodes = node_set[n.children[c2].label]
					for i in c1_nodes:
						if i.label not in dist_matrix:
							dist_matrix[i.label] = {}
						if i.label not in LCAs:
							LCAs[i.label] = {}
						for j in c2_nodes:
							if j.label not in dist_matrix:
								dist_matrix[j.label] = {}
							if j.label not in LCAs:
								LCAs[j.label] = {}
							dist_matrix[i.label][j.label] = height[i.label] + height[j.label]
							dist_matrix[j.label][i.label] = height[i.label] + height[j.label]
							LCAs[j.label][i.label] = n.label
							LCAs[i.label][j.label] = n.label
			continue

		if n.label not in dist_matrix:
			dist_matrix[n.label] = {}

		if n.label not in LCAs:
			LCAs[n.label] = {}

		node_set[n.label] = [n]
		for c1 in range(len(n.children) -1):
			c1_nodes = node_set[n.children[c1].label]
			for c2 in range(c1+1,len(n.children)):
				c2_nodes = node_set[n.children[c2].label]
				for i in c1_nodes:
					if i.label not in dist_matrix:
						dist_matrix[i.label] = {}

					if i.label not in LCAs:
						LCAs[i.label] = {}

					dist_matrix[i.label][n.label] = height[i.label]
					dist_matrix[n.label][i.label] = height[i.label]

					LCAs[i.label][n.label] = n.label
					LCAs[n.label][i.label] = n.label

					for j in c2_nodes:
						if j.label not in dist_matrix:
							dist_matrix[j.label] = {}

						if j.label not in LCAs:
							LCAs[j.label] = {}

						dist_matrix[j.label][n.label] = height[j.label]
						dist_matrix[n.label][j.label] = height[j.label]

						LCAs[j.label][n.label] = n.label
						LCAs[n.label][j.label] = n.label

						dist_matrix[i.label][j.label] = height[i.label] + height[j.label]
						dist_matrix[j.label][i.label] = height[i.label] + height[j.label]

						LCAs[j.label][i.label] = n.label
						LCAs[i.label][j.label] = n.label

				for i in c1_nodes:
					height[i.label] += n.edge_length
					node_set[n.label] += [i]
				for j in c2_nodes:
					height[j.label] += n.edge_length
					node_set[n.label] += [j]
				height[n.label] =  n.edge_length

	return dist_matrix, LCAs


def perm_test(dist_matrix, k, rest_anchors, original_stat, num = 1000):
	all_prob = []
	for _ in range(num):
		last_anchor = np.random.choice([d for d in dist_matrix if d not in rest_anchors])
		# anchors = np.random.choice([d for d in dist_matrix], size=k, replace=False)
		p = np.random.rand(k)
		p /= np.sum(p)

		min_p = k-1
		# last_anchor = anchors[-1]

		# rest_anchors = [a for a in anchors if a != last_anchor]
		min_dist = min([dist_matrix[int(last_anchor)][a] for a in rest_anchors])

		choices = set()
		for a in rest_anchors:
			choices.update([j for j in dist_matrix[a] if dist_matrix[a][j] <= min_dist])

		prob = len([c for c in choices if c not in rest_anchors]) / (len(dist_matrix) - k + 1)
		prob2 = 1 - (1 - p[min_p])**(k-1)
		# min_prob = min(min_prob, prob * prob2)

		all_prob.append(prob * prob2)
	all_prob = np.array(all_prob)
	return len(all_prob[all_prob < original_stat]) / num

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

def find_individual_y(p, x, ymean, anchors, tree_obj):
	pass


def save_jplace(all_rounds, tree_obj, file):
	result = {}
	tree_str, label_dict = extended_newick(tree_obj)
	result["metadata"] = {"invocation": " ".join(sys.argv),
						"software": f"DecoDiPhy {__version__}",
						"repository" : "https://github.com/shayesteh99/DecoDiPhy"}

	result["fields"] = ["edge_num", "abundance", "x", "y"]

	f = all_rounds[-1]
	placements = []
	for i in range(f['k']):
		p = {}
		p['n'] = "q"+str(i+1)
		index = label_dict[f['anchors'][i]]
		
		p['p'] = [[index, f['p'][i], f['x'][i], f['y']]]
		placements.append(p)

	result['placements'] = placements
	result["tree"] = tree_str

	with open(file, 'w') as f:
		f.write(json.dumps(result, sort_keys=True, indent=4))
		f.write("\n")
		f.close()
	# print(result)


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--tree', required=True, help="Input tree")
	parser.add_argument('-d', '--distances', required=True, help="Distance file")
	parser.add_argument('-s', '--seed', required=False, default=1142, help="Random Seed")
	parser.add_argument('-f', '--fix_k', required=False, default='0', help="Fix k")
	parser.add_argument('-k', '--k', required=False, help="Number of query taxa if fix_k==1 or Number of max queries if fix_k==0")
	parser.add_argument('-q', '--quick', required=False, default='1', help="Quick?")
	parser.add_argument('-r', '--radius', required=False, default=2, help="Radius")
	# parser.add_argument('-p', '--pval', required=False, default=0.1, help="P-value threshold")
	parser.add_argument('-m', '--method', required=False, choices=['hill', 'exhaustive', 'closest', 'closest_iterative', 'reverse_hill'], default='hill', help="Search method")
	parser.add_argument('-o', '--outdir', required=False, default='', help="Output dir")
	parser.add_argument('--min_p', required=False, default=1e-2, help="Minimum abundance")
	parser.add_argument('--warm_start', required=False, help="Path to the input file")

	args = parser.parse_args()

	random.seed(a=int(args.seed))
	np.random.seed(int(args.seed))

	# tracemalloc.start()
	start = time.time()

	# pval_thresh = float(args.pval)

	with open(args.tree,'r') as f:
		tree = f.read().strip().split("\n")[0]
	tree_obj = read_tree_newick(tree)
	if not __label_tree__(tree_obj):
		print("Input tree was not labeled!")
		with open(os.path.join(args.outdir, "labeled_tree.trees"), 'w') as f:
			f.write(tree_obj.newick())


	with open(args.distances, 'r') as f:
		lines = f.readlines()
		lines = [l.split() for l in lines]
		distances = {l[0]:float(l[1]) for l in lines}

	pruned = preprocess_input(tree_obj, distances)
	d, l, C, D, index_to_node, node_to_index, index_to_leaf = get_input_matrices(tree_obj, distances)
	D_T = D.T
	C_T = C.T
	l_T = l.reshape(-1)

	#only for p_value
	# dist_matrix = get_dist_matrix_new(tree_obj, node_to_index)

	all_rounds = []
	# _, goal_obj, _, _, _ = reverse_hill_climbing(d, l, C, D, index_to_node)
	# print(goal_obj)
	# return

	print("created matrices")

	if args.fix_k == '1':
		if not args.k:
			raise ValueError("Please specify the number of queries k.")
		k = int(args.k)
		k = min(k, len(index_to_node))

	if args.method != 'hill':
		if args.method == "exhaustive":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = exhaustive_search(d, l, C, D, index_to_node, k)
		elif args.method == "closest":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = k_closest_leaves(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k)
		elif args.method == "closest_iterative":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = k_closest_leaves_iterative(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k)
		elif args.method == "reverse_hill":
			opt_anchors, opt_obj, opt_p, opt_x, opt_y = reverse_hill_climbing_fixed_k(d, l, C, D, index_to_node, k, min_p = float(args.min_p))
		# elif args.method == "hill":
		# 	if args.quick == '1':
		# 		node_neighbors = get_neighbors(tree_obj, node_to_index, int(args.radius))
		# 	opt_anchors = None
		# 	for i in range(1, k+1):
		# 		opt_anchors, opt_obj, opt_p, opt_x, opt_y, last_anchor = hill_climbing(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors, quick = args.quick, neighbors = node_neighbors)
		print("opt_anchors: ", opt_anchors)
		print("opt_x: ", opt_x)
		print("opt_p: ", opt_p)
		print("opt_y: ", opt_y)

		round_info = {}
		round_info["k"] = len(opt_anchors)
		round_info["rounds"] = "final"
		round_info["loss"] = float(opt_obj)
		round_info["anchors"] = [index_to_node[i] for i in opt_anchors]
		round_info["p"] = list(opt_p)
		round_info["x"] = list(opt_x)
		round_info["y"] = float(opt_y)
		all_rounds.append(round_info)


	elif args.method == 'hill':
		opt_anchors = None
		prev_obj = None
		all_obj = []
		all_x = []
		all_p = []
		all_y = []
		all_anchors = []

		p_thresh = float(args.min_p)

		node_neighbors = None
		if args.quick == '1':
			node_neighbors = get_neighbors(tree_obj, node_to_index, int(args.radius))
		elif args.quick == '2':
			node_neighbors = get_immediate_neighbors(tree_obj, node_to_index)

		start_k = 1
		if args.warm_start:
			with open(args.warm_start, "r") as f:
				all_rounds = json.load(f)
			opt_anchors = np.array([node_to_index[a] for a in all_rounds[-1]['anchors']])
			start_k = len(opt_anchors) + 1

		end_k = len(index_to_node)
		if args.k:
			k = int(args.k)
			end_k = k+1

		if start_k >= end_k:
			if args.warm_start:
				with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
					json.dump(all_rounds, f)

			raise ValueError("k too small!")
			

		for i in range(start_k, end_k):
			print("+" * 200)
			print(i)

			#build cvxpy problem
			# problem = build_cvxpy_problem(D.shape[0], i)
			if opt_anchors is None:
				opt_anchors = np.random.choice([i for i in range(len(index_to_node))], i, replace = False)
			else:
				new = np.random.choice([i for i in range(len(index_to_node)) if i not in opt_anchors], i-len(opt_anchors), replace = False)
				opt_anchors = np.array(list(opt_anchors) + list(new))

			solver = AnchorOSQP(d, l_T, C_T, D_T, opt_anchors)

			# if args.parallel == '1':
			# 	opt_anchors, opt_obj, opt_p, opt_x, opt_y, round_info = hill_climbing_concurrent(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors, quick = args.quick, neighbors = node_neighbors)
			# else:
			# if args.quick == '2':
			# 	opt_anchors, opt_obj, opt_p, opt_x, opt_y, round_info = hill_climbing_even_faster(d, l, C, D, index_to_node, i, all_rounds, initial_anchors = opt_anchors, neighbors = node_neighbors, max_radius = int(args.radius))
			# else: 
			opt_anchors, opt_obj, opt_p, opt_x, opt_y, round_info = hill_climbing(d, l_T, C_T, D_T, index_to_node, i, all_rounds, solver, anchors = opt_anchors, quick = args.quick, neighbors = node_neighbors)
			print("opt_x: ", opt_x)
			print("opt_p: ", opt_p)
			print("opt_y: ", opt_y)

			# # if min(opt_p) < 1/len(index_to_node) and i > 2:
			# p_thresh = -np.log(1-0.25) / (i**2)
			# if i <= 2:
			# 	p_thresh = 1e-3 
			if min(opt_p) < p_thresh and args.fix_k == '0':
				print("stopped by p0")
				break



			# p_val = 1
			# if i >= 2:
			# 	min_prob = 1
			# 	print([index_to_node[a] for a in opt_anchors], opt_p)
			# 	for j in range(len(opt_anchors)):
			# 		last_anchor = opt_anchors[j]
			# 		min_p = j

			# 		rest_anchors = [a for a in opt_anchors if a != last_anchor]

			# 		min_dist = min([dist_matrix[last_anchor][a] for a in rest_anchors])

			# 		choices = set()
			# 		for a in rest_anchors:
			# 			choices.update([j for j in dist_matrix[a] if dist_matrix[a][j] <= min_dist])

			# 		prob = len([c for c in choices if c not in rest_anchors]) / (len(dist_matrix) - i + 1)
			# 		# prob2 = 1 - np.exp(-(i**2) * opt_p[min_p])
			# 		prob2 = 1 - (1 - opt_p[min_p])**(i-1)
			# 		min_prob = min(min_prob, prob*prob2)

			# 	p_val = perm_test(dist_matrix, i, rest_anchors, min_prob)
			# 	print("p val: ", p_val)

			# 	if p_val < pval_thresh and args.fix_k == '0':
			# 		print("stopped by p_val")
			# 		round_info = all_rounds[-1]
			# 		round_info['p_value'] = p_val
			# 		break

			# round_info = all_rounds[-1]
			# round_info['p_value'] = p_val

			with open(os.path.join(args.outdir, "all_rounds_" + str(i) + ".json"), 'w') as f:
				json.dump(all_rounds, f)
			if os.path.exists(os.path.join(args.outdir, "all_rounds_" + str(i-1) + ".json")):
				os.remove(os.path.join(args.outdir, "all_rounds_" + str(i-1) + ".json"))

			all_obj.append(opt_obj)
			all_x.append(opt_x)
			all_p.append(opt_p)
			all_y.append(opt_y)
			all_anchors.append(opt_anchors)
			if opt_obj < 1e-10 and args.fix_k == '0':
				print("stopped by obj0")
				break
			# if prev_obj and (prev_obj - opt_obj)/prev_obj < 0.05:
			# 	print("stopped by obj_imp")
			# 	break

			prev_obj = opt_obj

		if len(all_anchors) == 0:
			end = time.time()

			round_info = {}
			round_info["k"] = 0
			round_info["rounds"] = "final"
			round_info["loss"] = "NA"
			round_info["anchors"] = []
			round_info["p"] = []
			round_info["x"] = []
			round_info["y"] = 0
			round_info["runtime"] = end - start
			all_rounds.append(round_info)

			with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
				json.dump(all_rounds, f)

			return

		opt_anchors = all_anchors[-1]
		opt_x = all_x[-1]
		opt_y = all_y[-1]
		opt_p = all_p[-1]

	if max(opt_x) > 1 - 1e-3:
		indices = np.where(opt_x > 1 - 1e-3)[0]
		for i in indices:
			label = index_to_node[opt_anchors[i]]
			node = tree_obj.label_to_node(selection='all')[label]
			if not node.parent.is_root():
				parent = node_to_index[node.parent.label]
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
					new_obj, new_p, new_x, new_y, _ = get_optimal_obj_osqp(d, l_T, C_T, D_T, new_anchors, index_to_node, len(new_anchors))
					# new_obj, new_p, new_x, new_y, _ = get_optimal_obj(d, l, C, D, new_anchors, index_to_node, len(new_anchors))

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
			node = tree_obj.label_to_node(selection='all')[label]
			for c in node.child_nodes():
				s = node_to_index[c.label]
				if s in opt_anchors:
					continue
				new_anchors = np.array(opt_anchors)
				new_anchors[i] = s
				# new_obj, new_p, new_x, new_y, _ = get_optimal_obj(d, l, C, D, new_anchors, index_to_node, len(new_anchors))
				new_obj, new_p, new_x, new_y, _ = get_optimal_obj_osqp(d, l_T, C_T, D_T, new_anchors, index_to_node, len(new_anchors))
				if new_x[i] > 1e-4 and np.fabs(new_x[i] - 1) > 1e-3:
					opt_anchors = new_anchors
					opt_x = new_x
					opt_y = new_y
					opt_p = new_p
					break

	root_edges = [node_to_index[c.label] for c in tree_obj.root.child_nodes()]
	if len(root_edges) == 2 and root_edges[0] in opt_anchors and root_edges[1] in opt_anchors:
		if args.fix_k == '0':
			list(opt_anchors).remove(root_edges[1])
			opt_anchors = np.array(opt_anchors)
			# opt_obj, opt_p, opt_x, opt_y, _ = get_optimal_obj(d, l, C, D, opt_anchors, index_to_node, len(opt_anchors))
			opt_obj, opt_p, opt_x, opt_y, _ = get_optimal_obj_osqp(d, l_T, C_T, D_T, opt_anchors, index_to_node, len(new_anchors))

	end = time.time()

	sorted_indices = np.argsort(opt_anchors)
	opt_anchors = opt_anchors[sorted_indices]
	opt_anchors = [index_to_node[i] for i in opt_anchors]
	opt_p = opt_p[sorted_indices]
	opt_x = opt_x[sorted_indices]

	if pruned:
		opt_anchors, opt_x = post_process_output(tree_obj, read_tree_newick(tree), opt_anchors, opt_x)


	round_info = {}
	round_info["k"] = len(opt_anchors)
	round_info["rounds"] = "final"
	round_info["loss"] = all_rounds[-1]["loss"]
	round_info["anchors"] = list(opt_anchors)
	round_info["p"] = list(opt_p)
	round_info["x"] = list(opt_x)
	round_info["y"] = float(opt_y)
	round_info["runtime"] = end - start
	all_rounds.append(round_info)


	print("optimal anchors: ", opt_anchors)
	print("optimal p: ", opt_p)
	print("optimal x: ", opt_x)
	print("optimal y: ", opt_y)

	print("Optimization Runtime: ", end - start) 

	# current, peak = tracemalloc.get_traced_memory()
	# print(f"Current memory: {current / 1e6:.2f} MB")
	# print(f"Peak memory: {peak / 1e6:.2f} MB")

	# tracemalloc.stop()

	with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
		json.dump(all_rounds, f)
	if os.path.exists(os.path.join(args.outdir, "all_rounds_" + str(len(opt_anchors)) + ".json")):
		os.remove(os.path.join(args.outdir, "all_rounds_" + str(len(opt_anchors)) + ".json"))

	save_jplace(all_rounds, tree_obj, os.path.join(args.outdir, "output.jplace"))


if __name__ == "__main__":
	main()  