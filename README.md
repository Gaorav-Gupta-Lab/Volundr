# Völundr

## Python package for CRISPR Based Synthetic Lethal Screens in Mammalian Cells

## Getting Started

These instructions should allow a functional copy of this package to be installed on your local computer.  It is highly recomended
that you read the user guide before attempting to use this program.

## Prerequisites

Völundr has been designed to run in a 64 bit Linux environment.  It has been tested on Scientific Linux, RHEL, and CentOS.  It will
not run on any Windows OS.  Völundr does the target searching step in parallel.  The more CPUs or threads available the faster the 
job will complete.

```
Minimum System Requirements:
    4 CPUs or threads
    20 Gb RAM
    ~5x the size of the compressed FASTQ file for disk space.

Required Python Libraries:
    Python ≥v3.4
    numpy
    scipy
    python-levenshtein
    python-magic
    natsort
    pathos
    wxPython (Needed if using GUI API)
```

## Installation

The quickest way to install is to clone or download this repository into a location that accessable to the user.  

Install using setup.py

```
Exact steps coming soon
```

## Usage
The recomended method to run is using a bash shell as described in the user guide.

Völundr can be run from the command line.  If installed by cloning or downloading the repository:
```
python3 /path/to/Volundr.py --options_file /path/to/optionsfile
```

If installed using setup.py:
```
python3 Volundr --options_file /path/to/optionsfile
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Authors

* **Dennis Simpson** - *Initial work* 

## Cite

[![PubMed](img/2318832.png)](https://www.ncbi.nlm.nih.gov/pubmed/28649649)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
