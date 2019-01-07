# Flexible Structural Alignment
Flexible protein structural alignment with a hierarchical segmentation technique.

One of the recurring problems in structural alignment is to compare structures
belonging to the same *fold* family but in different conformations. The
classical algorithms are inappropriate because they consider both proteins to align as
rigid bodies. This problem is also critical when the proteins belong to
the same folding family but were obtained by experimental methods in
different conformations.

This project aims to achieve a new method of structural alignment by pair
of two proteins. The novelty brought by this method is to cut out one of the two
proteins into small subunits beforehand, in order to allow superimposition by "flexible" alignment and facilitate the alignment of proteins in different conformations.

## Installation

### Clone the repository
```
git clone https://github.com/gabrielctn/flexible_structural_alignment.git
cd flexible_structural_alignment
```

### Requirements

1. A linux distribution.

2. Install the few required **python packages** :

```
pip3 install -r requirements.txt

# This command will install the following modules:
# biopython == 1.72
# schema == 0.6.8
# docopt==0.6.2
```


## Run the program

`flex` takes in input two protein structures as two PDB files. The proteins should have only one chain A. In any case the program reindexes both files to make them have 1 chain A starting with residue index 1.


### Toy example

In the terminal run:
```
./flex data/1MBA.pdb data/101M.pdb
```

#### Get help

```
$./flex -h

Usage:
    ./flex PDB_FILE_1 PDB_FILE_2

Arguments:
    PDB_FILE_1         Path to the first PDB file to align
    PDB_FILE_2         Path to the second PDB file to align

Options:
    -h, --help        Show this
```

## Author

- [Gabriel Cretin](https://github.com/gabrielctn)

## License

This project is licensed under the MIT License.
