# contig-puller
Extracts contigs harbouring a target gene / sequence from assemblies and aligns the contigs around the target sequence

## Author
Jason Kwong (@kwongjc)  
GitHub: [kwongj](https://github.com/kwongj)  

## Dependencies
* [Python 2.x](https://www.python.org/downloads/)
* [BioPython](http://biopython.org/wiki/Main_Page)
* [Prokka](https://github.com/tseemann/prokka)

## Usage
`$ contig-puller.py -h`  
```
usage: 
  contig-puller.py --db DBFILE --out OUTFILE --anno ANNOFILE [OPTIONS] FASTA-1 FASTA-2 ... FASTA-N

Extracts contigs containing a target sequence (eg. gene) from several multi-FASTA assemblies 
 and aligns the contigs at the target sequence

positional arguments:
  FASTA            FASTA assemblies to search

optional arguments:
  -h, --help       show this help message and exit
  --db DBFILE      target sequence (FASTA) to search for (REQUIRED to create
                   BLAST database)
  --id ID          percentage identity cutoff (default=100)
  --out OUTFILE    output file (default=contigs.gbk)
  --anno ANNOFILE  reference proteins.faa file to annotate from (optional |
                   requires Prokka to be installed)
  --cpus CPUS      number of cpus to use (default=8)
  --version        show program's version number and exit
```

**Requires:**
* Target sequence in FASTA format to search for eg. blaKPC-2
* Assemblies to search for target sequence

**Options:**
* Set BLAST search %identity cutoff (default = 100) using `--id`
* Specify proteins.faa file to annotate from using `--anno`  
(recommended to visualise surrounding context of target sequence)

## Bugs
Please submit via the GitHub issues page: [https://github.com/kwongj/contig-puller/issues](https://github.com/kwongj/contig-puller/issues)  

## Software Licence
GPLv2: [https://github.com/kwongj/contig-puller/blob/master/LICENSE](https://github.com/kwongj/contig-puller/blob/master/LICENSE)

## Other links
* [Prokka](https://github.com/tseemann/prokka)
