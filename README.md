# Differential Splicing Statistics

The repository contains two scripts:
- functions_ds.R: A R-script containing all functions necessary for the evaluation.
- plots_ds.Rmd: A R-notebook which calls theses functions and is used for visualization.

Both files should be located in the same directory to work properly.

## Usage
functions_ds.R provides the following functions:
For reading the groundtruth:
- read_info(path): Takes a path to the info file generated by the AS-Simulator and reads it.
- read_annotation(path): Takes a path to the annotation file generated by the AS-Simulator and reads it.
- load_introns(read_annotation_output,TxDb): Takes the object returned by read_annotation and the genome annotatoin (as TxDb object) and creates additional event location annotations.

For processing the different tool outputs:
- read_ds(path,tool): Takes a path to the output directory of a tool and the tool name to read its output.
- pre_filter_type(read_ds_output,read_info_output): Takes a parsed tool output and the read info, to remove events and gene with unknown types.
- count_ds(read_ds_output,tool,...): Takes a parsed tool output (potentially filtered by pre_filter_type) to count tp, fp, fn for three different levels, which get returned as a data.table.

Supported tool names for the functions read_ds and count_ds are:
cash, eventpointer, aspli, majiq, spladder, edger, psisigma and junctionseq
The names are used to call the right parser internally

## Plotting
plots_ds.Rmd can be used to visualize the results.

