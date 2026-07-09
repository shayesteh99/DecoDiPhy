import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import argparse
import time
import json
import random
import warnings

import numpy as np
from treeswift import read_tree_newick

warnings.filterwarnings("ignore", category=RuntimeWarning)

from .__init__ import __version__
from .tree_utils import label_tree, preprocess_input, post_process_output, check_identifiability, compute_individual_ys
from .matrix import get_input_matrices
from .neighbors import get_neighbors, get_immediate_neighbors
from .search import hill_climbing, exhaustive_search, k_closest_leaves, k_closest_leaves_iterative, reverse_hill_climbing_fixed_k
from .output import save_jplace
from .optimize import AnchorOSQP, AnchorCVXPY, solve_osqp_once
from .read_jplace import jplace_to_distance, assignments_to_distance


def main():
    print("DecoDiPhy", __version__)

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, help="Input tree")
    parser.add_argument('-d', '--distances', required=False, help="Distance file")
    parser.add_argument('-a', '--assignments', required=False, help="A placement file with two columns: labels and counts")
    parser.add_argument('-j', '--jplace', required=False, help="jplace file")
    parser.add_argument('-s', '--seed', required=False, default=1142, help="Random Seed")
    parser.add_argument('-f', '--fix_k', required=False, default='0', help="Fix k")
    parser.add_argument('-k', '--k', required=False, help="Number of query taxa if fix_k==1 or Number of max queries if fix_k==0")
    parser.add_argument('-q', '--quick', required=False, default='1', help="Quick?")
    parser.add_argument('-r', '--radius', required=False, default=2, help="Radius")
    parser.add_argument('-m', '--method', required=False, choices=['hill', 'exhaustive', 'closest', 'closest_iterative', 'reverse_hill'], default='hill', help="Search method")
    parser.add_argument('-o', '--outdir', required=False, default='', help="Output dir")
    parser.add_argument('--min_p', required=False, help="Minimum abundance.")
    parser.add_argument('--warm_start', required=False, help="Path to the input file")
    parser.add_argument('--score', required=False, help="Score a set of placements")
    parser.add_argument('--optimizer', required=False, default='osqp', choices=['osqp', 'cvxpy'])

    args = parser.parse_args()

    random.seed(a=int(args.seed))
    np.random.seed(int(args.seed))

    if args.optimizer == "cvxpy":
        SolverClass = AnchorCVXPY
    elif args.optimizer == "osqp":
        SolverClass = AnchorOSQP
    else:
        raise ValueError(f"Unknown optimizer: {args.optimizer}")

    start = time.time()

    if args.min_p:
        p_thresh = float(args.min_p)

    if args.jplace:
        tree_obj, distances, read_count = jplace_to_distance(args.jplace)
        print("Total number of reads:", read_count)

        if not args.min_p:
            p_thresh = 1000 / read_count

        print("Input tree was not labeled! The labeled tree is saved in ", os.path.join(args.outdir, "labeled_tree.trees"))
        with open(os.path.join(args.outdir, "labeled_tree.trees"), 'w') as f:
            f.write(tree_obj.newick())

    elif args.assignments:
        if not args.tree:
            raise ValueError("You need to provide an input tree with your assignments.")

        with open(args.tree, 'r') as f:
            tree = f.read().strip().split("\n")[0]
        tree_obj = read_tree_newick(tree)
        if not label_tree(tree_obj):
            print("Input tree was not labeled!")
            with open(os.path.join(args.outdir, "labeled_tree.trees"), 'w') as f:
                f.write(tree_obj.newick())

        with open(args.assignments, 'r') as f:
            lines = f.readlines()
            lines = [l.split() for l in lines]
            assignments = {l[0]: float(l[1]) for l in lines}

        distances, read_count = assignments_to_distance(assignments, tree_obj)
        print("Total number of reads:", read_count)

        if not args.min_p and read_count >= 1000:
            p_thresh = 1000 / read_count
        else:
            p_thresh = 0.01

    else:
        if not args.tree or not args.distances:
            raise ValueError("You need either a jplace file or a tree and the average distances as input.")

        with open(args.tree, 'r') as f:
            tree = f.read().strip().split("\n")[0]
        tree_obj = read_tree_newick(tree)
        if not label_tree(tree_obj):
            print("Input tree was not labeled!")
            with open(os.path.join(args.outdir, "labeled_tree.trees"), 'w') as f:
                f.write(tree_obj.newick())

        with open(args.distances, 'r') as f:
            lines = f.readlines()
            lines = [l.split() for l in lines]
            distances = {l[0]: float(l[1]) for l in lines}

        if not args.min_p:
            p_thresh = 0.01

    pruned = preprocess_input(tree_obj, distances)
    d, l, C, D, index_to_node, node_to_index, index_to_leaf, branch_lengths, terminal_lengths = get_input_matrices(tree_obj, distances)
    D_T = D.T
    C_T = C.T
    l_T = l.reshape(-1)

    if args.score:
        with open(args.score, 'r') as f:
            lines = f.readlines()
            anchors = [l.split()[0] for l in lines]
            is_known = check_identifiability(tree_obj, anchors)
            print(is_known)
            solver = AnchorCVXPY(d, l_T, C_T, D_T, [node_to_index[a] for a in anchors])
            new_obj, new_p, new_x, new_y, _, _ = solver.solve([node_to_index[a] for a in anchors])

            for i in range(len(anchors)):
                print(anchors[i], "\t", new_p[i], "\t", new_x[i])

            print("loss: ", new_obj, is_known)
            print("anchors: ", anchors)
            print("optimal p: ", new_p)
            print("optimal x: ", new_x)
            print("optimal y: ", new_y)

            end = time.time()
            print("Optimization Runtime: ", end - start)
            return

    all_rounds = []

    print("Input matrices created!")

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
            opt_anchors, opt_obj, opt_p, opt_x, opt_y = reverse_hill_climbing_fixed_k(d, l, C, D, index_to_node, k, min_p=p_thresh)

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
            end_k = k + 1

        if start_k >= end_k:
            if args.warm_start:
                with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
                    json.dump(all_rounds, f)
            raise ValueError("k too small!")

        for i in range(start_k, end_k):
            print("=" * 200)
            print("k = ", i)

            if opt_anchors is None:
                opt_anchors = np.random.choice([i for i in range(len(index_to_node))], i, replace=False)
            else:
                new = np.random.choice([i for i in range(len(index_to_node)) if i not in opt_anchors], i - len(opt_anchors), replace=False)
                opt_anchors = np.array(list(opt_anchors) + list(new))

            solver = SolverClass(d, l_T, C_T, D_T, opt_anchors)

            opt_anchors, opt_obj, opt_p, opt_x, opt_y, round_info = hill_climbing(
                d, l_T, C_T, D_T, index_to_node, i, solver,
                anchors=opt_anchors, quick=args.quick, neighbors=node_neighbors,
            )

            opt_y = compute_individual_ys([index_to_node[i] for i in opt_anchors], opt_x, opt_y, terminal_lengths, branch_lengths)
            round_info['y'] = list(opt_y)

            all_rounds.append(round_info)
            print("Optimal Placements: ", *[index_to_node[i] for i in opt_anchors], sep = " ")
            print("Loss: ", opt_obj)
            # print("opt_x: ", opt_x)
            # print("opt_p: ", opt_p)
            # print("opt_y: ", opt_y)

            if min(opt_p) < p_thresh and args.fix_k == '0':
                print("Stopped by low abundance!")
                break

            with open(os.path.join(args.outdir, "all_rounds_" + str(i) + ".json"), 'w') as f:
                json.dump(all_rounds, f)
            if os.path.exists(os.path.join(args.outdir, "all_rounds_" + str(i - 1) + ".json")):
                os.remove(os.path.join(args.outdir, "all_rounds_" + str(i - 1) + ".json"))

            all_obj.append(opt_obj)
            all_x.append(opt_x)
            all_p.append(opt_p)
            all_y.append(opt_y)
            all_anchors.append(opt_anchors)
            if opt_obj < 1e-10 and args.fix_k == '0':
                print("Stopped by near-zero residue!")
                break

            prev_obj = opt_obj

        if len(all_anchors) == 0:
            end = time.time()
            round_info = {
                "k": 0,
                "rounds": "final",
                "loss": "NA",
                "anchors": [],
                "p": [],
                "x": [],
                "y": 0,
                "runtime": end - start,
            }
            all_rounds.append(round_info)
            with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
                json.dump(all_rounds, f, indent=4)
            return

        opt_anchors = all_anchors[-1]
        opt_x = all_x[-1]
        opt_y = all_y[-1]
        opt_p = all_p[-1]
        opt_obj = all_obj[-1]

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
                    new_obj, new_p, new_x, new_y, _ = solve_osqp_once(d, l_T, C_T, D_T, new_anchors, index_to_node, len(new_anchors))
                    new_y = compute_individual_ys([index_to_node[i] for i in new_anchors], new_x, new_y, terminal_lengths, branch_lengths)
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
                new_obj, new_p, new_x, new_y, _ = solve_osqp_once(d, l_T, C_T, D_T, new_anchors, index_to_node, len(new_anchors))
                new_y = compute_individual_ys([index_to_node[i] for i in new_anchors], new_x, new_y, terminal_lengths, branch_lengths)
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
            opt_obj, opt_p, opt_x, opt_y, _ = solve_osqp_once(d, l_T, C_T, D_T, opt_anchors, index_to_node, len(new_anchors))
            opt_y = compute_individual_ys([index_to_node[i] for i in opt_anchors], opt_x, opt_y, terminal_lengths, branch_lengths)

    end = time.time()

    sorted_indices = np.argsort(opt_p)[::-1]
    opt_anchors = opt_anchors[sorted_indices]
    opt_anchors = [index_to_node[i] for i in opt_anchors]
    opt_p = opt_p[sorted_indices]
    opt_x = opt_x[sorted_indices]

    if pruned:
        opt_anchors, opt_x = post_process_output(tree_obj, read_tree_newick(tree), opt_anchors, opt_x)

    # opt_y = compute_individual_ys(opt_anchors, opt_x, opt_y, terminal_lengths, branch_lengths)
    # print(opt_y)
    round_info = {
        "k": len(opt_anchors),
        "rounds": "final",
        "loss": all_rounds[-1]["loss"],
        "anchors": list(opt_anchors),
        "p": list(opt_p),
        "x": list(opt_x),
        "y": list(opt_y),
        "runtime": end - start,
    }
    all_rounds.append(round_info)

    print("=" * 200)
    print("Optimal Placments: ")
    for i in range(len(opt_anchors)):
        print("Edge: ", opt_anchors[i], ", abundance: ", opt_p[i], ", distal length: ", opt_x[i], ", terminal length: ", opt_y[i])
    # print("Optimal Placments: ", opt_anchors)
    # print("optimal p: ", opt_p)
    # print("optimal x: ", opt_x)
    # print("optimal y: ", opt_y)
    print("Optimization Runtime: ", end - start)

    with open(os.path.join(args.outdir, "all_rounds.json"), 'w') as f:
        json.dump(all_rounds, f, indent=4)
    if os.path.exists(os.path.join(args.outdir, "all_rounds_" + str(len(opt_anchors)) + ".json")):
        os.remove(os.path.join(args.outdir, "all_rounds_" + str(len(opt_anchors)) + ".json"))

    save_jplace(all_rounds, tree_obj, os.path.join(args.outdir, "output.jplace"))
