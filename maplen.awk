#!/usr/bin/awk -f

{
  chr=$1
  if ($4 > 0) {
    if (min[chr] == "") {
      min[chr] = $3
    }
    max[chr] = $3
  }
}

END {
  for (chr in min) {
    nchr += 1
    chr_length = max[chr] - min[chr]
    if (min_chr == "" || chr_length < min_chr) {
      min_chr = chr_length
    }
    total_length += chr_length
  }
  if (total_length < 100) {
    print "Total map length is < 100: rescaling to convert to cM\n"
    total_length *= 100
    min_chr *= 100
  }
  printf "Total_length\tMin_chrom_length\tNum_chroms\n"
  printf "%.4f\t%.4f\t\t\t%d\n", total_length, min_chr, nchr
}
