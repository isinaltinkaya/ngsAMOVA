# ngsAMOVA/plot
# Functions for plotting ngsAMOVA results

## Usage

```R
source("plot/ngsAMOVA_plot.R")
```

## `plot_tree` - plot tree from newick file

```R
plot_tree("amovaput.newick")
```

## `get_tree_from_newick()` - get tree object from newick file

```R
tree<-get_tree_from_newick("amovaput.newick")
```

You can then use the tree object to plot the tree using the `plot` function in the `ape` package.

```R
plot(tree)
```

___


# Misc
## Get block size parameter from ngsLD output

```sh
Rscript get_block_size.R --ld_files ld_file_list --fit_level 3 --seed 42 --decay_threshold 99
```
