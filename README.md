# Project Title

IBIS is a fast IBD Segment calling algorithm aimed at large, unphased genome datasets.

## Getting Started

Downloading this repository will give you the compiled c++ executable "ibis" and this is runnable by itself, though it may not be compiled optimally for all environments. If you want to compile it yourself, you will need to use the included Makefile.

### Prerequisites

Download the Makefile and the included IBIS.cc file, and run:

```
make
```
in the same directory. You will need the genetio library from another williamslab repository:
https://github.com/williamslab/genetio

## Supported input formats

IBIS requires PLINK .bed, .bim, and .fam format data to run. PLINK can be run with --make-bed to convert many other forms of genetic data into this file format.

##IBIS Usage

IBIS accepts its input .bed, .bim, and .fam files in one of two ways:

* First Three Arguments: [bed file] [bim file] [fam file]         
	*Specifies the plink format files for the data by specific name. Must be first 3 arguments.
or
* -b [prefix] or -bfile [prefix]         
	*Specifies the prefix to be used with prefix.bed, prefix.bim, and prefix.fam for the plink format input.
        *Does not need to be first argument
IBIS options:
*-mL or -min_l <value>            
	*Specify minimum length for acceptible segments to output.
	*Defaults to 7 centimorgans.
*-er or -errorRate <value>        
	*specify acceptible error rate in a segment before considering it false.                           
	*Defaults to .004 errors per marker
*-f or file <filename>           
	*Specify output file.
	*Defaults to ibis<thread number>.seg and will output a separate output file for each thread.
*-m or -merge                    
	*allows internal merging of segments over gaps caused by potential errors.
*-2 or -ibd2                     
	*enable ibd2 analyses
*-mt <value>                     
	*set an error threshold for leniency toward local spikes in error rate
*-threads <value>                
	*set the number of threads available to IBIS for parallel processing.


## License

This project is licensed under the GPL-3.0 - see the [LICENSE.md](LICENSE.md) file for details
