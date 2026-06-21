import sys
import json

from .__init__ import __version__
from .jutil import extended_newick


def save_jplace(all_rounds, tree_obj, file):
    result = {}
    tree_str, label_dict = extended_newick(tree_obj)
    result["metadata"] = {
        "invocation": " ".join(sys.argv),
        "software": "DecoDiPhy",
        "version": __version__,
        "repository": "https://github.com/shayesteh99/DecoDiPhy",
    }
    result["fields"] = ["edge_num", "abundance", "x", "y"]

    f = all_rounds[-1]
    placements = []
    for i in range(f['k']):
        p = {}
        p['n'] = "q" + str(i + 1)
        index = label_dict[f['anchors'][i]]
        p['p'] = [[index, f['p'][i], f['x'][i], f['y'][i]]]
        placements.append(p)

    result['placements'] = placements
    result["tree"] = tree_str

    with open(file, 'w') as f:
        f.write(json.dumps(result, sort_keys=True, indent=4))
        f.write("\n")
        f.close()
