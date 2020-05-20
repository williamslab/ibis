# IBIS

IBIS is a fast IBD Segment calling algorithm aimed at large, unphased genetic datasets.

### Compiling IBIS

First, ensure zlib, including developmental headers, is installed (e.g., the `zlib1g-dev` package on Ubuntu).

Next, clone the repository by running

    git clone --recurse-submodules https://github.com/williamslab/ibis.git

(Alternatively, `git clone [repo]` followed by `git submodule update --init` in the cloned directory will do the same as the above.)

Now, compile by running

    make

in the repository directory (i.e., run `cd ibis` then `make`).

To pull IBIS updates, use

    git pull
    git submodule update --remote

### Steps for running IBIS

1. Convert input data into PLINK binary format (.bed, .bim, and .fam).
  * Running [PLINK](https://www.cog-genomics.org/plink2/) with `--make-bed` enables conversion from many other forms of genetic data into this file format.

2. Insert a genetic map into the bim file using [add-map-plink.pl](add-map-plink.pl).
  * A good human genetic map is the [HapMap II map](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/) ([http link](http://bit.ly/33s6XQG)). As of this writing, the latest version for build 37 is available [here](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz) ([http link](http://bit.ly/2QuvNu9)).
  * Example `add-map-plink.pl` command using bash:
```
./add-map-plink.pl my.bim [map directory]/genetic_map_GRCh37_chr{1..22}.txt > new.bim
```
  * this uses `my.bim` to create a file `new.bim` with the genetic map inserted. (For non-bash environments, supply the file for each chromosome after the bim file: `./add-map-plink.pl my.bim genetic_map_GRCh37_chr1.txt genetic_map_GRCh37_chr2.txt genetic_map_GRCh37_chr3.txt ...`.)
  * NOT RECOMMENDED: If you do not add a genetic map to the bim file, and you have 0 for all genetic positions in the input, IBIS will use a genetic map based on the physical positions in the input, treating 1Mb as 1cM.

3. Run IBIS using the specifications described below.

### IBIS Usage

IBIS accepts its input .bed, .bim, and .fam files in one of two ways:

* First Three Arguments: `[bed file] [bim file] [fam file]`
	* Specifies the plink format files for the data by specific name. Must be first 3 arguments.
Example:
```
./ibis test1-chr1.bed test1-chr1.bim test1-chr1.fam -min_l 7 -mt 500 -er .004 -f test1Out
```
or
* `-b [prefix]` or `-bfile [prefix]`
	* Specifies that files with name prefix.bed, prefix.bim, and prefix.fam should be read. Does not need to be first argument.

```
./ibis -bfile test1-chr1 -min_l 7 -mt 500 -er .004 -f test1Out
```
### IBIS Options:

#### Execution options:
* `-2` or `-ibd2`
	* Enable IBD2 analyses
* `-hbd`
	* Enable HBD segment detection
* `-chr <value>`
	* Specify a single chromosome for IBIS to process when given an input file with multiple chromosomes.
* `-t <value>` or `-threads <value>`
	* Set the number of threads available to IBIS for parallel processing.
* `-noConvert`
	* Prevent IBIS from attempting to convert putative Morgan genetic positions to centiMorgans.
	* IBIS makes this conversion if any input chromosome is <= 6 genetic units in length, `-noConvert` disables.
* `-maxDist <value>`
	* Set a maximum separation distance between SNPs in the input map. This is described in more detail in the IBIS paper.
	* Defaults to being inactive.

#### IBD threshold parameters:
* `-er <value>` or `-errorRate <value>`
	* specify acceptable error rate in a segment before considering it false.
	* Defaults to .004 errors per marker
* `-mL <value>` or `-min_l <value>`
	* Specify minimum length for acceptable segments to output.
	* Defaults to 7 centiMorgans.
* `-mt <value>`
	* Specify a minimum number of markers required for a segment to be printed.
	* Defaults to 448 markers

#### IBD2 threshold parameters: (use with -2 or -ibd2)
* `-er2 <value>` or `-errorRate2 <value>`
	* specify acceptible error rate in an IBD2 segment before considering it false.
	* Defaults to .008 errors per marker
* `-mL2 <value>` or `-min_l2 <value>`
	* Specify minimum length for acceptable IBD2 segments to output.
	* Defaults to 2 centiMorgans.
* `-mt2 <value>`
	* Specify a minimum number of markers required for an IBD2 segment to be printed.
	* Defaults to 192 markers

#### HBD threshold parameters: (use with -hbd)
* `-erH <value>` or `-errorRateH <value>`
	* specify acceptable error rate in an HBD segment before considering it false.
	* Defaults to .008 errors per marker
* `-mLH <value>` or `-min_lH <value>`
	* Specify minimum length for acceptable HBD segments to output.
	* Defaults to 3 centiMorgans.
* `-mtH <value>`
	* Specify a minimum number of markers required for an HBD segment to be printed.
	* Defaults to 192 markers

#### Output controls:
* `-f <filename>` or `-o <filename>` or `-file <filename>`
	* Specify output file prefix.
	* Defaults to `ibis`, resulting in ibis.seg, ibis.coef (if using `-printCoef`), ibis.hbd (with `-hbd`), and ibis.incoef (with `-hbd` and `-printCoef`).
* `-bin` or `-binary`
	* Have the program print the segment data file in binary format. Use bseg2seg to interpret.
* `-gzip`
	* Have the program output gzipped segment and (optionally) coef files.
* `-noFamID`
	* Have IBIS use only the individual ID in its output notations.
	* Defaults to `<fam ID>:<indiv ID>`

#### Kinship and inbreeding coefficient file options:
* `-printCoef`
	* Have the program print a .coef file and a .incoef (with -hbd) file.
* `-a <value>`
	* Set a different supplemental kinship coefficient factor to add for each pair.
	* Defaults to 0.00138.
* `-d <value>` or `-degree <value>`
	* Set a minimum degree of relatedness. IBD segments and coefficients for pairs with more distant degrees will not be printed to the output.
	* Defaults to including all degrees.
	* Mutually exclusive with `-c`
* `-c <value>`
	* Set a minimum kinship coefficient. IBD segments and coefficients for pairs with lower kinships will not be printed to the output.
	* Defaults to 0.
	* Mutually exclusive with `-d`

### IBIS Output

IBIS produces a .seg or .bseg file and, when using `-printCoef`, a .coef file. Additionally, when using `-hbd`, it prints a .hbd file with homozygous by descent segments, and with the addition of `-printCoef`, a .incoef file.

Segment file format:
```
sample1 sample2 chrom phys_start_pos phys_end_pos IBD_type genetic_start_pos genetic_end_pos genetic_seg_length marker_count error_count error_density
```
* `IBD_type` can be either IBD1 or IBD2
* `error_count` and `error_density` are negative numbers for IBD1 segments that precede IBD2 segments. The error information in them is not specifically tracked by IBIS.

If `-bin` is employed, the .bseg output will not be human readable, and bseg2seg can be used to convert .bseg files to .seg format.

Any number of .bseg files can be provided for conversion.
Example:
```
./bseg2seg test1-chr1.fam test1Out.chrom1.bseg test1Out.chrom2.bseg test1Out.chrom3.bseg ...
```


Coef file format:
```
sample1 sample2 kinship_coefficient IBD2_fraction segment_count degree_of_relatedness
```
* Unrelated individuals have a degree of -1
* The supplemental coefficient factor (`-a`) is included in the given kinship coefficients and for determining the degrees.
* When -2 is not used, IBD1 and IBD2 are not distinguished from each other, and the threshold used for degree classification is based on IBD0 fraction. In this case, the supplemental coefficient factor is added to the IBD0 fraction.

HBD file format:
```
sample_id chrom phys_start_pos phys_end_pos HBD_type genetic_start_pos genetic_end_pos genetic_seg_length marker_count error_count error_density
```
* `HBD_type` is simply `HBD`

Incoef file format:
```
sample_id inbreeding_coefficient segment_count
```

### seg2coef

The `-printCoef` option requires IBIS to analyze genome-wide data, but it is possible to produce a coef file from a set of .seg or .bseg files. This would be useful when a user runs IBIS on each chromosome independently.

`seg2coef` takes the following options:

    ./seg2coef [total map length cM] [fam file] [seg/bseg files ...]

The total map length in cM is available in the .bim file. The maplen.awk script calculates this in the following way:

    ./maplen.awk [bim files ...]

**Example seg2coef execution** using data subdivided into chromosomes as `data[chr].{bed,bim,fam}`, with IBIS segments in `output[chr].seg`:

    ./maplen.awk data{1..22}.bim
    ./seg2coef [Total_length] data1.fam output{1..22}.seg > output.coef

the `[Total_length]` argument should be the value output from the `maplen.awk` command. This produces `output.coef`.

## License

This project is licensed under the GPL-3.0 - see the [LICENSE](LICENSE) file for details
