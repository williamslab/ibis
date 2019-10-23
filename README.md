# IBIS

IBIS is a fast IBD Segment calling algorithm aimed at large, unphased genome datasets.

### Compiling IBIS

Download the Makefile and the included IBIS.cc file, and run:

```
make
```
in the same directory. You will need the genetio library from another williamslab repository:
https://github.com/williamslab/genetio

IBIS also requires the zlib and libpthread libraries.
### Supported input formats

IBIS requires PLINK .bed, .bim, and .fam format data to run. PLINK can be run with --make-bed to convert many other forms of genetic data into this file format.

### IBIS Usage

IBIS is designed to process one chromosome at a time. This single chromosome can be specified with the -chr option, or alternatively input files can be provided with only one chromosome.

IBIS accepts its input .bed, .bim, and .fam files in one of two ways:

* First Three Arguments: [bed file] [bim file] [fam file]         
	* Specifies the plink format files for the data by specific name. Must be first 3 arguments.
Example:
```
./ibis test1-chr1.bed test1-chr1.bim test1-chr1.fam -min_l 7 -mt 500 -er .004 -f test1Out
```
or
* -b [prefix] or -bfile [prefix]         
	* Specifies the prefix to be used with prefix.bed, prefix.bim, and prefix.fam for the plink format input.
	* Does not need to be first argument

```
./ibis -bfile test1-chr1 -min_l 7 -mt 500 -er .004 -f test1Out
```
IBIS options:

* -chr \<value\>
	* Specify a single chromosome for IBIS to process when given an input file with multiple chromosomes.
	* IBIS will fail if this is not specified and there are multiple chromosomes in the input files.
* -mL or -min_l \<value\>            
	* Specify minimum length for acceptible segments to output.
	* Defaults to 7 centimorgans.
* -er or -errorRate \<value\>        
	* specify acceptible error rate in a segment before considering it false.                           
	* Defaults to .004 errors per marker
* -f or file \<filename\>           
	* Specify output file.
	* Defaults to ibis\<thread number\>.seg and will output a separate output file for each thread.
* -m or -merge                    
	* allows internal merging of segments over gaps caused by potential errors.
* -2 or -ibd2                     
	* enable ibd2 analyses
* -mt \<value\>                     
	* set an error threshold for leniency toward local spikes in error rate
* -threads \<value\>                
	* set the number of threads available to IBIS for parallel processing.

## License

This project is licensed under the GPL-3.0 - see the [LICENSE](LICENSE) file for details
