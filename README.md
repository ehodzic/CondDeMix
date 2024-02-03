# CondDeMix
Condition-dependent deconvolution of bulk expression data

### System Requirements
- make (version 3.81 or higher)
- g++, supporting c++11 or later (GCC version 5.2 or higher is used for linking with GUROBI in the Makefile)
- Gurobi Optimizer

### Compiling `CondDeMix`
In the `Makefile`, set `GUROBIROOT` to the path of your root Gurobi folder.

Simply run `make` command in the root CondDeMix folder. It will create the executables.

### Running `CondDeMix`
**Usage:**
```sh
./CondDeMix -m bulkfile -b referencefile -s sigfile -c sigName
```

| Parameters | Description |
| ------ | ------ |
| `-m` | The bulk expression file, in form of a tab-separated data matrix (genes x samples). Data is organized so that the first row is the header row. The header row contains a name for the first column, which is followed by all the sample names. Every subsequent row begins with the name of a gene, then followed by the expression of that gene in all the samples listed in the header. |
| `-b` | The single-cell reference file, in form of a tab-separated data matrix (genes x cell types). Data is organized so that the first row is the header row. The header row contains a name for the first column, which is followed by all the cell-type names. Every subsequent row begins with the name of a gene, then followed by the counts of that gene in all the cell types listed in the header. |
| `-s` | The mutational signature file, in form of a tab-separated data matrix (samples x signature) with two columns. Data is organized so that the first row is the header row. The header row contains a name for the first column, which is followed by the mutational signature name. Every subsequent row begins with the name of a sample, then followed by the signature strengths in that sample for the mutational signatures listed in the header row. |
| `-c` (optional) | A string indicating the name of the signature column (matching the name in the header of the signature matrix -s) to be used in the model. If this parameter is not provided, and the matrix in `-s` has only one column, then that column's name will be used. If the `-s` matrix has multiple columns and this parameter is not provided, an error will be reported and the program will not run. |
| `-k` (optional) | Name of the file containing a list of marker genes, each on a separate line. A header row is assumed to exist, with the column names, which is ignored when reading this file. Each subsequent line is assumed to contain tab-separated entries denoting, in this exact order: the cell type name (string), gene name (string), distance (real number), and the number of cell types that the marker is expressed in (integer). |
| `-f` (optional) | The path of the output folder, in which all the output files are placed. If this parameter is not provided, it will default to `./output/testRun`. |

#### Example

Command to run:
```sh
./CondDeMix -m data/bulkExpression.tsv -b data/referenceMatrix.tsv -k data/markerGenes.tsv -s data/signatureStrengths.tsv -c "SBS4" -f output/exampleRun -v 1
```
