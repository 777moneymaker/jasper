# JASPER
![JASPER LOGO](https://github.com/777moneymaker/jasper/blob/main/logo.png)
JASPER is a free bioinformatics tools for predicting virus hosts. 

------------

Jasper uses a bunch of bioinformatics tools to prediction virus hosts. It includes genome-genome alignment, CRISPR spacers analyzation, tRNA analyzation and more.


## Requirements
### Naming convention
**Jasper** depends on good file naming convention. The best is to use sequence ID as file name, e.x. `NC_008876.fna`. Software will use this id to name every temp file that needs to be created and also it will use this ID in results file.

If you put multiple contigs in a single file, there is no problem with that. Just be sure that every contig is in it's right file. **Jasper** repairs every file it read, by default naming it `<id from filename>|<#contig>` e.x.:
```
>NC_000856|1
ATGCT....
>NC_000856|2
ATGCA....
# and so on
```
So even if you have, for instance, one genome in your file, then **Jasper** will change it's id to `<id from filename>|1`.

------------

### Python Libraries
```
Biopython
Numpy
Pandas
Pytz
Six
```
### Additional software
```
NCBI-Blast+
PILER-CR
tRNAscan-SE
```
## Installation
##### *I will make a simple installation script in my free time.*

------------


**JASPER** uses additional software. It calls every program with `subprocess` so every program that is stated in above should be installed and added to `$PATH`. 

* To install NCBI-Blast+ use `sudo apt install ncbi-blast+`
* To install PILER-CR go [here](http://www.drive5.com/pilercr/), dowload compiled software, move somwhere and add to `$PATH`.
* To install tRNAscan-SE use `sudo apt install trnascan-se`

By defaults, tRNAscan-SE install as `trnascan-1.4`. Be sure to use `trnascan` alias in your bash. To do it:
```
vim ~/.bashrc
    alias trnascan="trnascan-1.4"
source ~/.bashrc
```

After that go to JASPER's main directory and:
```
python3 setup.py install
```

## Usage
JASPER uses bunch of arguments. A lot of parameters are BLAST parameters and can be configured with JSON file and passed to JASPER.

##### Basic usage
```
jasper --virus path/to/virus/dir --host /path/to/host/dir
```

##### Advanced usage
```
TODO
```


## Additional information
[PILER-CR](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-18)
[NCBI-BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

## License
[MIT](https://choosealicense.com/licenses/mit/)
