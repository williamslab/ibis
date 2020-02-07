# IBIS

IBIS is a fast IBD Segment calling algorithm aimed at large, unphased genome datasets.

### Compiling IBIS

Download the Makefile and the included IBIS.cc file, and run:

```
make
```
in the same directory. You will need the genetio library from:
https://github.com/williamslab/genetio


Include the genetio.a file's location in your LD_LIBRARY_PATH environment variable.

IBIS also requires zlib and OpenMP.
### Supported input formats

IBIS requires PLINK .bed, .bim, and .fam format data to run. PLINK can be run with --make-bed to convert many other forms of genetic data into this file format.

### IBIS Usage

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
### IBIS Options:

* -chr \<value\>
	* Specify a single chromosome for IBIS to process when given an input file with multiple chromosomes.
* -mL or -min_l \<value\>            
	* Specify minimum length for acceptible segments to output.
	* Defaults to 7 centimorgans.
* -mL2 or -min_l2 \<value\>
	* Specify minimum length for acceptible IBD2 segments to output.
	* Defaults to 2 centimorgans.
* -er or -errorRate \<value\>        
	* specify acceptible error rate in a segment before considering it false.                           
	* Defaults to .004 errors per marker
* -er2 or -errorRate2 \<value\>
	*specify acceptible error rate in an IBD2 segment before considering it false.
* -f or file \<filename\>           
	* Specify output file prefix.
	* Defaults to ibis, resulting in ibis.seg.\<thread number\> and ibis.coef.\<thread number\> and will output a separate output file for each thread.
* -2 or -ibd2                     
	* enable ibd2 analyses
* -mt \<value\>                     
	* Specify a minimum number of markers required for a segment to be printed.
* -mt2 \<value\> 
	* Specify a minimum number of markers required for an IBD2 segment to be printed.
* -threads \<value\>                
	* set the number of threads available to IBIS for parallel processing.
* -gzip
	* Have the progrem output gzipped segment and coef files.
* -printCoef
	* Have the progrem print the .coef files. 
* -c \<value\>
	* Set a minimum kinship coefficient. Pairs with lower kinship will not be printed to the output.
* -d \<value\>
	* Set a minimum degree of relatedness. Pairs with more distant degrees will not be printed to the output.
* -a \<value\>
	* Set a different supplemental coefficient factor to add to pairs.
	* Defaults to 0.00138
* -noFamID
	* Have IBIS use only the individual ID in its output notations.
	* Defaults to \<fam ID\>:\<indiv ID\>
* -force
	* Prevent IBIS from attempting to convert putative Morgan genetic positions to centiMorgans by multiplying these by 100
	* IBIS makes this conversion if any input chromosome is <= 6 genetic units in length, -force disables
### IBIS Output

Ibis produces 2 files for each thread: a .seg file and a .coef file.

Thread file format:
```
sample1 sample2 chrom phys_start_pos phys_end_pos IBD_type gen_start_pos gen_end_pos seg_length marker_count error_count error_density
```
* IBD_type can be either IBD1 or IBD2
* error_count and error_density are negative numbers for IBD1 segments that are being printed due to preceeding printed IBD2 segments, and the error information in them is not available to IBIS.


Coef file format:
```
sample1 sample2 kinship_coefficient IBD2_fraction segment_Count degree_of_relationship
```
* Unrelated individuals have a given degree of -1 
* The supplemental coefficient factor is included in the given kinship coefficients and degrees.
* When -2 is not used, and IBD1 and IBD2 are not distinguished from each other, the threshold used for degree classification is based on IBD0 fraction, and the supplemental coefficient factor is added to the IBD0 fraction by multiplying it by 4 to convert it into genome fraction units.


## License

This project is licensed under the GPL-3.0 - see the [LICENSE](LICENSE) file for details
