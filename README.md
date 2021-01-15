# JASPER
![JASPER LOGO](https://github.com/777moneymaker/jasper/blob/main/logo.png)

------------

JASPER is a free bioinformatics tools for predicting virus hosts. 
JASPER uses a bunch of bioinformatics tools to prediction virus hosts. It includes genome-genome alignment, CRISPR spacers analyzation, tRNA analyzation and more.
JASPER contains few, independent modules `blast`, `crispr`, `trna - not included now`, `merge`.

## Requirements
### Naming convention
**Jasper** depends on good file naming convention. The best is to use sequence ID as file name, e.x. `NC_008876.fna`. Software will use this id to name every temp file that needs to be created and also it will use this ID in results file.

**WARNING** It's not the best idea to use `|` char in your filename and also in sequence header. Just use normal fasta naming like `>NC_00876 additional_info more_additional_info`.

If you put multiple contigs in a single file, there is no problem with that. Just be sure that every contig is in it's right file. **Jasper** repairs every file it read, by default naming it `<id from filename>|<#contig>` e.x.:
```
>NC_000856|1
ATGCT....
>NC_000856|2
ATGCA....
# and so on
```
So even if you have, for instance, one genome in your file, then **Jasper** will change it's id to `<id from filename>|1`.

**EXTENSIONS!** Jasper uses input files that ends with (fa, fna, fasta) only!

------------

### Additional software
```
NCBI-Blast+
PILER-CR
tRNAscan-SE # not included in current release
```
## Installation
##### *I will make a simple installation script in my free time.*

------------

**JASPER** uses additional software. It calls every program with `subprocess` so every program that is stated in above should be installed and added to `$PATH`. 

On Ubuntu:
* To install NCBI-Blast+ use `sudo apt install ncbi-blast+`
* To install PILER-CR go [here](http://www.drive5.com/pilercr/), dowload compiled software, move somwhere and add to `$PATH` under name `pilercr`.
* To install tRNAscan-SE use `sudo apt install trnascan-se`

Source code for additional software:
* NCBI-Blast+: [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* PILER-CR [here](http://www.drive5.com/pilercr/)
* tRNAscan-SE [here](http://lowelab.ucsc.edu/tRNAscan-SE/)

**Remember to install everything and add it to path**

After that go to JASPER's main directory and:
```
python setup.py install
```
## Tests
If you want to test, go to proj directory and type `python -m unittest discover`.

## Usage
JASPER uses bunch of arguments. A lot of parameters are BLAST parameters and can be configured with JSON file and passed to JASPER.

##### Basic usage
```
python3 jasper[.py] blast --virus path/to/virus/dir --create-db host_db --host /path/to/host/dir
python3 jasper[.py] crispr --host path/to/host/dir --create-db vir_db --host /path/to/vir/dir
python3 jasper[.py] trna --args # soon
python3 jasper[.py] merge results1.csv results2.csv results3.csv [-w --weigths] 0.6 0.2 0.2 --output final_results.csv 
```

For more check `--help` on jasper individual modules: `python3 jasper[.py] {blast,crispr,trna,merge} --help`

##### Blast config
You can provide blast config as as a `*.json` file.
Every module uses different task so there are few arguments that are forbidden:
`[query', 'db', 'outfmt', 'max_target_seqs', 'num_alignments']`

## References
* Edgar, R.C. (2007) [PILER-CR: fast and accurate identification of CRISPR repeats](http://www.ncbi.nlm.nih.gov/pubmed/17239253), BMC Bioinformatics, Jan 20;8:18
* Fichant and Burks, J. Mol. Biol. (1991) Identification of tRNA genes in genomic DNA, 220:659-671.
* [NCBI-BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

## License
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)
