# SGM-Clustering-Programs
Tools for Evaluating Performance of SGM Clustering Alternatives
This is a framework for building your own set of benchmarks.

## Compilation
SGM and SGM_Loop can be compiled simply by using g++ assuming you have boost-libs (or its distribution-specific equivalent) installed and available in your PATH.

`g++ -O3 SGM.cpp -o SGM`

## Use

SGM requires a fasta file and a primary, secondary, and tertiary threshold.  An optional fifth parameter will assign a user-defined prefix to all output.

`SGM example-genome.fasta 0.90 0.95 0.95 example-output-prefix`

#### Thresholds

The primary threshold controls MJSD-based splicing.  By iterating over the entire genome, SGM will identify the maximum MJSD between two segments at all possible splice-points if one were to split the genome.  If the significance of this MJSD is at or above the first threshold, the splice is made.  The same process is repeated recursively for all resultant segments until no segment can be further segmented at or above the threshold.

The second threshold is the initial, adjacent clustering.  Assuming hyper-segmentation, segments are recursively scanned and recombined if the significance of their MJSD is below this threshold until all segments can be iterated upon without a recombination event.

The final threshold is only applicable to hierarchical and standard clustering methods of SGM.
In standard clustering, segments are iterated over sequentially, and recombination events are generated on a first-observed basis.  For example, the first segment is compared to the second through last segments until a recombination event occurs.  If recombination occurs, the first segment becomes a combination of these two segments, and the segment it was combined with is removed from the list, and the iteration process restarts.  If no recombination event occurs, the second segment is compared to the third through final segments, and so on.

Once all segments can be iterated over without any recombination, standard clustering stops.

For hierarchical clustering, all possible combinations of segments and their MJSD is calculated. 
Each round, the most significant MJSD that is also below the provided third threhsold is combined.  This is repeated, until no more combinations below the threshold can be combined, regardless of their distance from each other within the genome.

### Modification
SGM can produce a variety of outputs containing the cluster sequences, including:
- Edge-Weight Precursor Files for MCL Clustering (See abcmaker.py for more information)
- Dissimilarity Matrix Files
- Similarity Matrix Files
- Fasta Files

These may be turned on/off within the CPP files for SGM (SGM.cpp, SGM_Loop.cpp) by commenting/uncommenting the output lines, e.g. changing

`netfilePrinter(outfile4,stage3segments);` 

to 

`//netfilePrinter(outfile4,stage3segments);`

Will stop SGM from producing the precursor files necessary for producing the ABC files used by the MCL clustering method.

Hierarchical clustering is provided by the clusterbeyondX function.
Standard (non-hierarchical) clustering is provided by the clusterbeyond function.  Adjust these depending on your needs.

### SGM_Loop

This program will iterate through a range of supplied parameters for the third threshold (if applicable) for each type of clustering method.  This is very useful for benchmarking.  Simply provide a Ranges.txt file containing the desired thresholds on individual lines.  A Ranges.txt file containing, for example, 0.001 and 0.002 on individual lines and combined with the following command:

`./SGM_Loop <Input Genome Fasta> 0.90 0.95 <Output File Name>`

Will produce clustering combinations of the following first, second, and third threshold parameters.

- 0.90 0.95 0.001
- 0.90 0.95 0.002

Modifying the clustering methods and output for SGM_Loop.cpp is the same as SGM.cpp.

### abcmaker.py

The output from the netfilePrinter function (.network files) are provided as input to this program to generate MCL-compatible ABC (.abc) files.
See the MCL documentation at https://micans.org/mcl/ for more information.

Edge weights are proportional to the significance of the MJSD between all segment pairs.
To lower computational resources / generate sparse networks, a threhsold parameter can be set to filter output edges.

`python3 abcmaker.py input.network`
