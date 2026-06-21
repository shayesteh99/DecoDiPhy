import sys
import os
import argparse
import json
import re
import gzip

import numpy as np
from treeswift import read_tree_newick, Node

from .tree_utils import label_tree, get_terminal_branch_lengths


def preprocess_tree(tree_obj, epsilon=1e-4):
    for n in tree_obj.traverse_preorder():
        if not n.is_root():
            if n.edge_length <= 0:
                n.edge_length = epsilon


def read_tree(tree):
    new_tree = re.sub(r"\{.*?\}", "", tree)
    tree_obj = read_tree_newick(new_tree)
    is_labeled = label_tree(tree_obj)
    preprocess_tree(tree_obj)
    return tree_obj, is_labeled


def merge_placements(placements, labels, is_weighted=False, mode="count"):
    branch_count = {}
    for p in placements:
        if len(p['p']) == 0:
            continue
        length = len(p['p'])
        total = length
        if is_weighted:
            total = sum(p['p'][i][4] for i in range(length))
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

    total = sum(branch_count[k] for k in branch_count)
    if mode == "count":
        return branch_count
    branch_count = {k: branch_count[k] / total for k in branch_count}
    return branch_count


def place_queries(tree, branch_count):
    tree_obj = read_tree_newick(tree)
    terminal_lengths = get_terminal_branch_lengths(tree_obj)

    for q in branch_count:
        node = tree_obj.label_to_node(selection='all')[q]
        tlength = terminal_lengths[q]
        if node.is_root():
            query = Node(label='q' + str(q), edge_length=tlength)
            parent = Node(label='p' + str(q), edge_length=0)
            node.edge_length = 0
            parent.add_child(node)
            parent.add_child(query)
            tree_obj.root = parent
        else:
            query = Node(label='q' + str(q), edge_length=tlength)
            parent = Node(label='p' + str(q), edge_length=node.edge_length / 2)
            node.edge_length = node.edge_length / 2
            nparent = node.parent
            nparent.remove_child(node)
            nparent.add_child(parent)
            parent.add_child(node)
            parent.add_child(query)
    return tree_obj


def compute_avg_distances(tree_obj, branch_count):
    distance_matrix = tree_obj.distance_matrix(leaf_labels=True)
    avg_matrix = {}
    query_labels = ['q' + str(q) for q in branch_count]

    for q in branch_count:
        d = distance_matrix['q' + str(q)]
        for l in d:
            if l not in query_labels:
                if l not in avg_matrix:
                    avg_matrix[l] = d[l] * branch_count[q]
                else:
                    avg_matrix[l] += d[l] * branch_count[q]

    return avg_matrix


def match_trees(tree_obj, query_tree, branch_count):
    id_to_clade = {}
    num_to_id = {}

    for n in tree_obj.traverse_preorder():
        n_leaves = [l.label for l in n.traverse_leaves()]
        id_to_clade[n.label] = set(n_leaves)

    for n in query_tree.traverse_preorder():
        n_leaves = [l.label for l in n.traverse_leaves()]
        set_leaves = set(n_leaves)
        is_l = False
        for k in id_to_clade:
            if id_to_clade[k] == set_leaves:
                num_to_id[k] = n.label
                is_l = True
                break
        if not is_l:
            print(n.label, set_leaves)

    labels = [n.label for n in tree_obj.traverse_postorder()]
    ids_to_labels = {i: labels[i] for i in range(len(labels))}

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
                branch_count[c] = p / len(child_ids)
            else:
                branch_count[c] += p / len(child_ids)
        del branch_count[root_id]

    for b in branch_count:
        new_branch_count[num_to_id[ids_to_labels[b]]] = branch_count[b]

    return new_branch_count


def jplace_to_distance(file):
    file_format = file.split(".")[-1]

    if file_format == "gz":
        with gzip.open(file, 'rt') as f:
            placements = json.load(f)
    else:
        with open(file, 'r') as f:
            placements = json.load(f)

    tree = placements['tree']
    tree_obj, is_labeled = read_tree(tree)

    labels = [n.label for n in tree_obj.traverse_postorder()]
    branch_count = merge_placements(placements['placements'], labels)

    total = sum(branch_count[k] for k in branch_count)
    branch_count = {k: branch_count[k] / total for k in branch_count}

    tree_with_queries = place_queries(tree_obj, branch_count)
    avg_matrix = compute_avg_distances(tree_with_queries, branch_count)

    return tree_obj, avg_matrix, total

def assignments_to_distance(assignments, tree):
    tree_labels = [n.label for n in tree.traverse_preorder() if not n.is_root()]
    assignments = {i:assignments[i] for i in assignments if i in tree_labels}

    total = sum([assignments[k] for k in assignments])
    assignments = {k: assignments[k] / total for k in assignments}
    tree_with_queries = place_queries(tree, assignments)

    avg_matrix = compute_avg_distances(tree_with_queries, assignments)

    return avg_matrix, total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help="Input Placements")
    parser.add_argument('-w', '--weighted', required=False, default='0', help="For multiplacements only")
    parser.add_argument('-m', '--mode', required=False, choices=['count', 'abundance'], default='count', help="mode")
    parser.add_argument('-o', '--output', required=True, help="Output file")
    parser.add_argument('-t', '--tree', required=False, help="labelled tree file")

    args = parser.parse_args()

    file_format = args.input.split(".")[-1]

    if file_format == "gz":
        with gzip.open(args.input, 'rt') as f:
            placements = json.load(f)
    else:
        with open(args.input, 'r') as f:
            placements = json.load(f)

    tree = placements['tree']
    tree_obj, is_labeled = read_tree(tree)

    if not is_labeled:
        tree_obj.write_tree_newick(args.tree)

    is_weighted = args.weighted == '1'

    labels = [n.label for n in tree_obj.traverse_postorder()]
    branch_count = merge_placements(placements['placements'], labels, is_weighted, args.mode)

    with open(args.output, 'w') as f:
        for l in branch_count:
            f.write(l + "\t" + str(branch_count[l]) + "\n")


if __name__ == "__main__":
    main()
