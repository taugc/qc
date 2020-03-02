![](/images/taugc.png) 
# TauGC Bioinformatics - Quality Control 

**Under Contruction!**

## Installation using miniconda:

  1- [Install miniconda](https://conda.io/miniconda.html)

  2 - Create conda env using Python3:
  ```
  $ conda create -n qc python=3
  ```
  3- In/Out conda the environment:
  ```
  $ conda activate qc
  $ conda deactivate
  ```
  4- Install main packages using conda:
  ```
   conda install -c bioconda --file requirements.txt
  ```

  5- Running:

  ```
  $ python qc_taugc.py -h
  ```

  ```
  usage: qc_taugc.py [-h] [-r R] id {SE,PE} {0,1} f

  QC TauGC: used to make the quality control from NGS fastq files and generate
  an automatic report.

  positional arguments:
    id          Sample name
    {SE,PE}     Choose the type of sequencing: (SE): Single-end (PE): Pair-end.
    {0,1}       Choose experiment: (0): RNAseq (1): Others
    f           Forward fastq.gz file.

  optional arguments:
    -h, --help  show this help message and exit
    -r R        Reverse fastq.gz file.
  ```

  ## Code of conduct
  Please, check our code of conduct [here](/docs/CODE_OF_CONDUCT.md).