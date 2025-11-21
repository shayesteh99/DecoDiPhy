# DecoDiPhy

## Installation 
You can clone DecoDiPhy from github:

```
git clone https://github.com/shayesteh99/DecoDiPhy.git
cd DecoDiPhy
```

## Exucation
To use DecoDiPhy, you can simply run:
```
python DecoDiPhy.py -d [INPUT DISTANCE FILE] -t [REFERENCE TREE] -o [OUTPUT DIR] > [LOG FILE]
```

Example:
```
python DecoDiPhy.py -t Example/pruned_tree.trees -d Example/true_distances.txt -o ./Example
```

Below is the description of all parameters of DecoDiPhy:
>- `--tree` or `-t`: The reference tree in the newick format. See `Example/pruned_tree.trees`.
>- `--distances` or `-d`: The file with the sample distance to all taxa in the reference tree. See `Example/true_distances.txt` for formatting.
>- `--min_p`: This is the minimum abundance of individual placements that DecoDiPhy will detect. This is a way to control the resolution of consolidation or number of placements by DecoDiPhy. Higher means earlier stopping and less placements. Lower means higher resolution and more placements. The default value is `0.01`, however, we sugggest setting it to `1000/|X|` where `|X|` is the total number if reads in your sample.
>- `--k` or `-k`: This is the maximum number of placements that DecoDiPhy will find. This is a more explicit way of controlling how many placements you want from DecoDiPhy. By setting `-k [K]`, DecoDiPhy will stop at `[K]` placements if none of the other stopping criteria holds (e.g. `min_p`). If you want DecoDiphy to find *exactly* `[K]` placements, use `-k [K]` with `-f 1`.
>- `--outdir ` or `-o`: This is the output directory where DecoDiPhy will use. By the end of running DecoDiPHy, you'll see two files in this directory: `output.jplace` is the final placements in jplace format, and `all_rounds.json` contains all the placements DecoDiPhy finds throughout its search in JSON format. It will be useful if you're wondering what placements DecoDiPhy would output with lower resolution. It is also useful if you decide you want DecoDiPhy to continue its search and find mpre placements. You don't need to start over, as it can be time extensive. You can simply use the `--warm_start [OUTPUT DIR]/all_rounds.json`, and DecoDiPhy will start from where it ended last time :)

### Using krepp output as input to DecoDiPhy:
krepp is a recent tool developed by Ali Osman Berk Şapci and Siavash Mirarab (See [preprint](https://www.biorxiv.org/content/10.1101/2025.01.20.633730v1)). You can visit [here](https://github.com/bo1929/krepp/tree/master) to check it out and learn how to use it. krepp will output distances to each reference taxa and placements of individual reads on the backbone tree, and tha't why it pairs very well with DecoDiPhy. In our experience, individual krepp plcacements result in more accurate consolidation with DecoDiPhy. If you have the krepp output in jplace format like `Example/krepp_multiplacements.jplace`, you can simply run:

```
python compute_distances.py -i Example/krepp_multiplacements.jplace -t Example/pruned_tree.trees -o Example/krepp_distances.txt
```

And now you can run DecoDiPhy on top of krepp distances:
```
python DecoDiPhy.py -t Example/pruned_tree.trees -d Example/krepp_distances.txt -o ./Example
```
