# JASPER
 **J**UST **A** **S**IMPLE **PRE**DICTOR.

![pypi](https://img.shields.io/pypi/v/jasper-vh.svg?branch=master)
![JASPER LOGO](https://github.com/777moneymaker/jasper/blob/main/logo.png?raw=true)

* JASPER is a free bioinformatics tools for predicting virus hosts. 
* JASPER uses a bunch of bioinformatics tools to prediction virus hosts. It includes genome-genome alignment, CRISPR spacers analyzation, tRNA analyzation and more.
* JASPER contains few, independent modules `blast`, `crispr`, `trna`, `wish`, `mash`, `merge`.

# Requirements
### Python 3.7
You need `Python >= 3.7` to use JASPER.

### Naming convention
**Jasper** depends on good file naming convention. The best is to use sequence ID as file name, e.x. `NC_008876.fna`. Software will use this id to name every temp file that needs to be created and also it will use this ID in results file.

**WARNING** It's not the best idea to use `|` char in your filename and also in sequence header. Just use normal fasta naming like `>NC_00876 additional_info more_additional_info`.

If you put multiple contigs in a single file, there is no problem with that. Just be sure that every contig is in it's right file. **Jasper** repairs every file, by default naming it `<id from filename>|<#contig>` e.x.:
```
>NC_000856|1
ATGCT....
>NC_000856|2
ATGCA....
# and so on
```
So even if you have, for instance, one genome in your file, then **Jasper** will change it's id to `<id from filename>|1`.

### Extensions
Jasper uses input files that ends with `[fa, fna, fasta]` only!

### Additional software
```
NCBI-Blast+
PILER-CR
WIsH
Mash
tRNAscan-SE
```

# Installation
**JASPER** uses additional software. It calls every program with `subprocess` so every program that is stated in above should be installed and added to `PATH`.

On Ubuntu:
* To install NCBI-Blast+ use `sudo apt install ncbi-blast+`
* To install PILER-CR go [here](http://www.drive5.com/pilercr/), download compiled software, move somewhere and add to `$PATH` under name `pilercr`.
* To install tRNAscan-SE go [here](http://lowelab.ucsc.edu/tRNAscan-SE/), download, compile, move somewhere and add to `$PATH` under name `tRNAscan-SE`. Remember that tRNAscan-SE needs `Infernal` to work properly.
* To install WIsH go [here](https://github.com/soedinglab/WIsH), download, compile, move somewhere and add to `$PATH` under name `WIsH`.
* To install Mash go [here](https://github.com/marbl/Mash), download release, move somewhere and add to `$PATH` under name `mash`.

Source code for additional software:
* NCBI-Blast+: [here](https://www.ncbi.nlm.nih.gov/books/NBK279671/)
* PILER-CR [here](http://www.drive5.com/pilercr/)
* tRNAscan-SE [here](http://lowelab.ucsc.edu/tRNAscan-SE/)
* WIsH [here](https://github.com/soedinglab/WIsH)
* Mash [here](https://github.com/marbl/Mash)

**Remember to install everything and add it to path**

**You can also download the script:** 
```
install_dependencies.sh - Linux
install_dependencies_osx.sh - Mac OS
```
After that go to JASPER's main directory and:
```
python setup.py install
```

### PATH

By defaults some pip on linux drops scripts to `~/.local/bin`. Add it to your `$PATH` at the end.
`export PATH="$HOME/.local/bin:$PATH"`
Now you're done, and you can start using `jasper-vh`.

### Installation summary
You want to have:
* installed jasper python package using `setuptools` or `pip`.
* installed each tool and added to `PATH`

# Tests
If you want to test, go to proj directory and type `python -m unittest discover`.
It's recommended to do that, since it performs tool check (ensures that user has all dependencies and proper python version).

# Usage
JASPER uses a bunch of arguments. A lot of parameters are BLAST parameters and can be configured with JSON file and passed to JASPER.
It's also recommended using `jasper` in empty directory. This ensures, that none of the user's file will be overwritten or damaged. Just do `mkdir jasper_results && cd jasper_results` and you're good to go. 

### Basic usage
```
jasper-vh blast --virus path/to/virus/dir --create-db host_db --host /path/to/host/dir --clear
jasper-vh crispr --host path/to/host/dir --create-db vir_db --host /path/to/vir/dir --clear
jasper-vh trna --host path/to/host/dir --virus /path/to/vir/dir --clear
jasper-vh wish --host path/to/host/dir --virus /path/to/vir/dir --clear
jasper-vh mash --host path/to/host/dir --virus /path/to/vir/dir --clear
jasper-vh merge {blast,crispr,trna,mash,wish}.csv --output final_results.csv
```
For more check `--help` on jasper individual modules: `jasper-vh  {blast,crispr,trna,wish,mash,merge} --help`

# Output
`JASPER` produces output in a special format: 

The resulting file is a csv file with additional lines that start with `#` (for easy parsing).

A resulting file is grouped by virus genome. For each viral genome, there are number of `Score` columns (number of score columns are eq to number of tools used/resulting files merged).

Under each group there is a STD (standard deviation) of that column, which indicates a level of variation. It gives an idea of how much the proposed host are different between each other in a group.

Sample:

|Virus    |Host     |blastScore|crisprScore|mashScore|wishScore|
|---------|---------|----------|-----------|---------|---------|
|NC_024389|NC_008531|294.0     |NaN        |0.537086 |-1.3626  |
|NC_024389|NC_008531|294.0     |NaN        |0.537086 |-1.3626  |
|NC_024389|NC_017331|NaN       |NaN        |NaN      |-1.35805 |
|# Std    |         |0.0       |NaN        |0.0      |0.002    |
|NC_024391|NC_009641|31180.0   |NaN        |0.760474 |-1.33066 |
|NC_024391|NC_017349|17892.0   |NaN        |0.789652 |-1.32703 |
|NC_024391|NC_017349|17892.0   |NaN        |0.789652 |-1.32703 |
|# Std    |         |6264.023  |NaN        |0.014    |0.002    |
|NC_024392|NC_018586|270.0     |28.0       |0.484772 |-1.36183 |
|NC_024392|NC_021827|270.0     |NaN        |NaN      |-1.36077 |
|NC_024392|NC_021830|270.0     |NaN        |NaN      |-1.3607  |
|NC_024392|NC_021839|270.0     |NaN        |NaN      |-1.36207 |
|NC_024392|NC_021840|270.0     |NaN        |NaN      |-1.35926 |
|NC_024392|NC_017547|74.0      |37.0       |0.517779 |-1.36224 |
|NC_024392|NC_021823|246.0     |37.0       |0.537086 |-1.36111 |
|NC_024392|NC_021823|246.0     |37.0       |0.537086 |-1.36111 |
|NC_024392|NC_010001|NaN       |NaN        |NaN      |-1.34728 |
|# Std    |         |63.37     |3.897      |0.021    |0.004    |

`Std` is equal to `NaN` only when the whole column is equal to `NaN` which means that there were no results in given tool for given hosts. 

# Blast config
You can provide blast config as a `*.json` file.
Every module uses a different task so there are few arguments that are forbidden:
`['query', 'db', 'outfmt', 'max_target_seqs', 'num_alignments']`

# References
* Edgar, R.C. (2007) [*PILER-CR: fast and accurate identification of CRISPR repeats*](http://www.ncbi.nlm.nih.gov/pubmed/17239253), BMC Bioinformatics, Jan 20;8:18
* Fichant and Burks, J. Mol. Biol. (1991) *Identification of tRNA genes in genomic DNA*, 220:659-671.
* Clovis Galiez, Matthias Siebert et al. *WIsH: who is the host? Predicting prokaryotichosts from metagenomic phage contigs*
* Ondov, B.D., Treangen, T.J., Melsted, P. et al. [*Mash: fast genome and metagenome distance estimation using MinHash.*](https://doi.org/10.1186/s13059-016-0997-x) Genome Biol 17, 132 (2016).
* [NCBI-BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

# License
[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)
