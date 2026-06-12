import time
import numpy as np


def hill_climbing(d, l_T, C_T, D_T, index_to_node, k, all_rounds, solver, anchors=None, quick='0', neighbors=None):
    opt_times = []
    opt_memos = []
    n = len(index_to_node)

    print(anchors)
    original, p, x, y, t, m = solver.solve(anchors)
    opt_times.append(t)
    opt_memos.append(m)
    min_p, min_x, min_y = p, x, y
    min_obj = original
    rounds = 0
    while True:
        init = time.time()
        rounds += 1
        print("=" * 200)
        print(rounds)
        min_obj = original
        og_anchors = np.array(anchors)

        for ind in range(k):
            i = k - ind - 1
            min_anchor = anchors[i]

            choices = range(n)
            if ind > 0 and k > 2 and quick == '1':
                choices = neighbors[anchors[i]]

            for v in choices:
                if v in anchors:
                    continue
                anchors[i] = v
                new_obj, new_p, new_x, new_y, t, _ = solver.solve(anchors)
                opt_times.append(t)
                if min_obj > new_obj:
                    min_obj = new_obj
                    min_anchor = v
                    min_p = new_p
                    min_x = new_x
                    min_y = new_y

            anchors[i] = min_anchor
            solver.solve(anchors)

        print(min_obj)
        end_round = time.time()
        round_info = {
            "k": k,
            "rounds": rounds,
            "loss": min_obj,
            "anchors": [index_to_node[i] for i in anchors],
            "p": list(min_p),
            "x": list(min_x),
            "y": float(min_y),
            "runtime": end_round - init,
            "opttime": np.mean(opt_times),
            "memory": np.mean(opt_memos),
        }
        all_rounds.append(round_info)
        if set(og_anchors) == set(anchors):
            return anchors, min_obj, min_p, min_x, min_y, round_info
        original = min_obj
        og_anchors = anchors

    return anchors, min_obj, min_p, min_x, min_y, round_info


def exhaustive_search(d, l, C, D, index_to_node, k):
    raise NotImplementedError(
        "exhaustive search requires the legacy cvxpy solver which is no longer available. Use --method hill instead."
    )


def k_closest_leaves(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k):
    raise NotImplementedError(
        "closest search requires the legacy cvxpy solver which is no longer available. Use --method hill instead."
    )


def k_closest_leaves_iterative(d, l, C, D, index_to_node, node_to_index, index_to_leaf, k):
    raise NotImplementedError(
        "closest_iterative search requires the legacy cvxpy solver which is no longer available. Use --method hill instead."
    )


def reverse_hill_climbing_fixed_k(d, l, C, D, index_to_node, k, min_p=0.01):
    raise NotImplementedError(
        "reverse_hill search requires the legacy cvxpy solver which is no longer available. Use --method hill instead."
    )
