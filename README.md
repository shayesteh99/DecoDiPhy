# DecoDiPhy
DecoDiPhy is a tool for consolidating noisy and erroneous sequence read information into accurate phylogenetic placements. All you need to run DecoDiPhy is a distance vector of your sample to each of your reference taxa, and a phylogeny on the same reference taxa. DecoDiPhy will find a multi-placement on this phylogeny that explains your sample.

## Installation 
DecoDiPhy is available on `bioconda`. You can install the latest version using the following commands:
```
conda create -n decodiphy python=3.11 decodiphy
conda activate decodiphy
```
Alternatively, you can install it directly from GitHub, and follow the steps below to create a local version:

```
conda create -n decodiphy python=3.11
conda activate decodiphy
git clone https://github.com/shayesteh99/DecoDiPhy.git
cd DecoDiPhy
pip install .
```

If the build is successful, running `decodiphy --help` will show you all the parameters and their use for decodiphy.

## Interactive Tutorial

If you're new to DecoDiPhy, we recommend starting with our interactive Google Colab tutorial. The tutorial walks through the complete workflow, including installation, preparing input files, running DecoDiPhy on a real metagenomic sample, interpreting the output files, visualizing the results, and controlling the number of placements.

**Launch the tutorial here:** **[DecoDiPhy Tutorial](https://github.com/shayesteh99/DecoDiPhy/blob/main/DecoDiPhy_tutorial.ipynb)**

The notebook is fully self-contained and can be run directly in your browser without installing DecoDiPhy locally.

## Execution
DecoDiPhy accepts multiple input formats. The easiest way to run it is to use the jplace output from krepp. krepp is a recent tool developed by Ali Osman Berk Şapci and Siavash Mirarab (See [preprint](https://www.biorxiv.org/content/10.1101/2025.01.20.633730v1)). You can visit [here](https://github.com/bo1929/krepp/tree/master) to check it out and learn how to use it. krepp will output distances to each reference taxa and placements of individual reads on the backbone tree, and that's why it pairs very well with DecoDiPhy. In our experience, individual krepp placements are more accurate for consolidation with DecoDiPhy. If you have the krepp output in jplace format like `Example/krepp_multiplacement.jplace`, you can simply run:
```
decodiphy -j [KREPP OUTPUT] -o [OUTPUT DIR] > [LOG FILE]
```

Example:
```
decodiphy -j Example/krepp_multiplacement.jplace -o ./Example
```

Alternatively, DecoDiPhy version 1.5.0 also supports genome assignments. You can use a simple assignment file with two columns: genome/node label and read count/abundance, and a backbone tree with the same labels. See `Example/krepp_assignments.txt` for how to prepare this file. Then, you can run:

```
decodiphy -a [INPUT ASSIGNMENT FILE] -t [REFERENCE TREE] -o [OUTPUT DIR] > [LOG FILE]
```

Example:
```
decodiphy -a Example/krepp_assignments.txt -t Example/reference_tree.trees -o ./Example
```

Below is the description of all parameters of DecoDiPhy:
>- `--min_p`: This is the minimum abundance of individual placements that DecoDiPhy will detect. This is a way to control the resolution of consolidation or the number of placements by DecoDiPhy. Higher means earlier stopping and fewer placements. Lower means higher resolution and more placements. The default value is `0.01`, however, we suggest setting it to `1000/|X|` where `|X|` is the total number of reads in your sample.
>- `--k` or `-k`: This is the maximum number of placements that DecoDiPhy will find. This is a more explicit way of controlling how many placements you want from DecoDiPhy. By setting `-k [K]`, DecoDiPhy will stop at `[K]` placements if none of the other stopping criteria hold (e.g. `min_p`). If you want DecoDiphy to find *exactly* `[K]` placements, use `-k [K]` with `-f 1`.
>- `--outdir ` or `-o`: This is the output directory where DecoDiPhy will use. By the end of running DecoDiPhy, you'll see two files in this directory: `output.jplace` is the final placements in jplace format, and `all_rounds.json` contains all the placements DecoDiPhy finds throughout its search in JSON format. It will be useful if you're wondering what placements DecoDiPhy would output with lower resolution. It is also useful if you decide you want DecoDiPhy to continue its search and find more placements. You don't need to start over, as it can be time extensive. You can simply use the `--warm_start [OUTPUT DIR]/all_rounds.json`, and DecoDiPhy will start from where it ended last time :)

### Output format
DecoDiPhy outputs two files: `output.jplace` and `all_rounds.json`. Below you can find more details about each of the files:

>- `output.jplace`: This jplace file contains the information for each placement. It follows the conventional jplace format proposed by [Matsen et al. (2012)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009). An example of the jplace output of decodiphy can be found in `Example/output.jplace`. For each placement, we report four fields: 
```
"fields": [
        "edge_num",
        "abundance",
        "x",
        "y"
    ]
```
where `edge_num` is the label of the placement edge in the backbone tree (the backbone tree in jplace format is labeled in a post-order traversal incrementally starting from 0), `abundance` is the abundance or proportion assigned to that placement, `x` is the relative placement position on the placement branch, and `y` is the terminal length. Currently, DecoDiPhy only outputs the average terminal lengths weighted by the abundances, so `y` will be the same for all placements.

>- `all_rounds.json`: This file contains more information about each round of DecoDiPhy optimization. It includes a list of JSON objects for each round of optimization. Below is an example of a JSON object for one round of optimization:

```
{
	"k": 3, 
	"rounds": 1, 
	"loss": 0.00018 
	"anchors": ["I104", "G000426325", "I9"], 
	"p": [0.72, 0.14, 0.14], 
	"x": [0.63, 0.35, 0.12], 
	"y": 0.0983, 
	"runtime": 0.5652, 
	"opttime": 0.0023
}
```
where `k` is the number of placements, `rounds` is the number of optimization rounds for that `k`, `loss` is the final loss value for the opimization function, `anchors` is the list of size `k` of placement edges labelled by initial labels of the backbone tree (if the internal nodes in the input tree are not labeled, DecoDiPhy will label them and save the labeled tree in `[OUTPUT DIR]/labeled_tree.trees`), `p` and `x` are lists of size `k` denoting the abundances and the relative placement positions, `y` is the average terminal length, and `runtime` and `opttime` are the running time of the total optimization and the small-prolem, respectively. The last JSON object in this list is labeled with `rounds="final"` and corresponds to the set of placements in the jplace file.

