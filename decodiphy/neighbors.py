from collections import deque


def get_neighbors(tree, node_to_index, radius=2):
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
