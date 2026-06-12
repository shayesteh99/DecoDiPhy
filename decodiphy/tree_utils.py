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

    distances = {d: distances[d] for d in distances if d in leaves}

    if len(leaves) < total:
        return True
    return False


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
        terminal_lengths[n.label] = sum_lengths[n.label] / leaf_counts[n.label] + n.edge_length / 2

    return terminal_lengths
