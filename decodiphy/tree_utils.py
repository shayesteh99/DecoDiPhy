import warnings


def is_number(x):
    try:
        float(x)
        return True
    except (ValueError, TypeError):
        return False


def label_tree(tree_obj, index=0):
    is_labeled = True
    labels = set()
    for node in tree_obj.traverse_preorder():
        if node.is_leaf():
            continue
        if not node.label or node.label in labels or is_number(node.label):
            is_labeled = False
            node.label = 'I' + str(index)
            index += 1
        labels.add(node.label)
    return is_labeled, index


def preprocess_input(tree_obj, distances):
    all_leaves = {l.label for l in tree_obj.traverse_leaves()}
    keep = {d for d in distances if d in all_leaves}
    # treeswift warns for every labeled unifurcation it suppresses
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        pruned_tree = tree_obj.extract_tree_with(keep)
    return pruned_tree, len(keep) < len(all_leaves)


def post_process_output(tree_obj, og_tree, opt_anchors, opt_x):
    new_anchors = []
    new_x = []
    for i in range(len(opt_anchors)):
        anchor = opt_anchors[i]
        x = opt_x[i]
        node = tree_obj.label_to_node(selection="all")[anchor]
        og_node = og_tree.label_to_node(selection="all")[anchor]
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


def check_identifiability(tree, anchors):
    label_to_node = tree.label_to_node(selection='all')
    for a in anchors:
        node1 = label_to_node[a]
        if node1.parent.label in anchors or node1.parent.parent.label in anchors:
            return True
    return False


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
        if n.is_root():
            terminal_lengths[n.label] = sum_lengths[n.label] / leaf_counts[n.label]
        else:
            terminal_lengths[n.label] = sum_lengths[n.label] / leaf_counts[n.label] + n.edge_length / 2

    return terminal_lengths


def compute_individual_ys(anchors, all_x, est_y, terminal_lengths, branch_lengths):
    all_y = []
    total_y = 0
    for i in range(len(anchors)):
        anchor = anchors[i]
        total_y += terminal_lengths[anchor] + all_x[i] * branch_lengths[anchor]
    offset = est_y/total_y

    for i in range(len(anchors)):
        anchor = anchors[i]
        y = offset * (terminal_lengths[anchor] + all_x[i] * branch_lengths[anchor])
        all_y.append(y)
    return all_y
