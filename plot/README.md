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

