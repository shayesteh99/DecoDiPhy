import numpy as np
from .tree_utils import get_terminal_branch_lengths


def get_input_matrices(tree_obj, distances):
    leaves = [n for n in tree_obj.traverse_leaves()]
    index_to_leaf = [l.label for l in leaves]
    leaf_to_index = {index_to_leaf[i]: i for i in range(len(index_to_leaf))}

    index_to_node = [n for n in tree_obj.labels(leaves=True, internal=True) if n != tree_obj.root.label]
    node_to_index = {index_to_node[i]: i for i in range(len(index_to_node))}

    C = np.zeros((len(index_to_leaf), len(index_to_node)))
    for l in tree_obj.traverse_leaves():
        node = l
        while not node.is_root():
            C[leaf_to_index[l.label], node_to_index[node.label]] = 1
            node = node.parent
    C[C == 0] = -1

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
            D[i, j] = dist_matrix[node][leaf]

    l = np.array([label_to_node[i].edge_length for i in index_to_node])
    l = l[:, np.newaxis]

    branch_lengths = {label_to_node[i].label:label_to_node[i].edge_length for i in index_to_node}
    terminal_lengths = get_terminal_branch_lengths(tree_obj)

    d = np.zeros(len(index_to_leaf))
    for i in range(len(index_to_leaf)):
        leaf = index_to_leaf[i]
        d[i] = distances[leaf]

    return d, l, C, D, index_to_node, node_to_index, index_to_leaf, branch_lengths, terminal_lengths
