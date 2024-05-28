## Description

A simple script that takes a list of UniProt IDs, and pulls residue-level annotations
from UniProt using the ID mapping tool via the public API.

As of now, it produces a .csv file with the following info:
* UniProt ID
* Residue ID
* Residue features (e.g. 'Natural variant', 'Modified residue' etc.)
* Curated variants (from humsavar)
* Subcellular locations
* GO terms (all and categorised)

## Getting Started

The program takes as input a list of UniProt IDs and outputs a .csv file 
with the annotated residues.

### Dependencies

No dependencies for now - only a modern version of Python (>=3.10).

### Installing

No installation needed, just clone the repo and it should be ready to run.

### Executing program

* Given that you have a list of UniProt IDs, you can run with something like:
```
python uniprot_annotation.py -i ids.list -o output.csv 
```

## Help

For usage help run:
```
python uniprot_annotation.py -h
```

## Authors

[Ioannis Riziotis](mailto:ioannis.riziotis@crick.ac.uk)
