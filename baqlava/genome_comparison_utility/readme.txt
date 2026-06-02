This utility provides two scripts that allow you to compare the genome/plasmid filtering levels across different filters that can be applied by BAQLaVa v1.1.0.

The first script, make_filtering_plot.py, can be used to visualize whether genomes are or are not included in the VGB for BAQLaVa at each filtering level (no filtering, V1 filtering, default, and conservative).

The second sript, make_filtering_and_mapping_plot.py, can be used to visualize the same information as script 1, with the added information of what your specific BAQLaVa marker or protein mappings looked like. This will help to inderstand whether you are seeing heavy mapping on genomes you would trust at one filtering level, but may not at others.

Usage: 

For each script, you should supply the VGB of interest (formatted exactly as shown in the BAQLaVa output, e.g. VGB_0000001), whether you are interested in the filtering of markers ("marker") or proteins ("protein"), and whether you are interested in the filtering of MGX ("MGX") or MTX ("MTX") profiles. 
You may additionally specify a path to the data files in this directory if you have moved them elsewhere. Otherwise, they are assumed to be in the same location as the script itself. 

For the second script, you will also provide the location of your tempfiles of interest. Use a path immediately upstream of your files of interest, as the script will search recursively to find all tempfiles contained in subdirectories. 

python make_filtering_plot.py [VGB_###] [marker|protein] [MGX|MTX] [path/to/data (optional)] 

make_filtering_and_mapping_plot.py [VGB_###] [marker|protein] [MGX|MTX] [path/to/tempfiles] [path/to/data (optional)]

Output:

The scripts will save out a plot with the name format automatically applied based on the parameters supplied: 

Script 1: <VGB>_<Markers/Proteins>_<MGX/MTX>_plasmid_filtering.png

Script 2: <VGB>_<Markers/Proteins>_<MGX/MTX>_mapping.png
