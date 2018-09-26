
### This repo is still under construction. The code is working, but not totally bug free, and the documentation is sparse, but the code has plenty of comments.

SAM_Reader is a Python module for working with many SAM/BAM formatted files. It is built upon [pysam](https://pysam.readthedocs.io/en/latest/api.html) which is great if you are working with one SAM/BAM file. 

## Quick Start Cookbook
**How many reads mapped to each organism across all files?**
```
import sam_reader

my_files = sam_reader('~/pathto/results/')
my_files.hits(write_file='~/new_file.csv')
```
**How many reads mapped to each organism from each file?**
**How many reads mapped to each posistion in an organism?**
**What does each base in the genome look like across all reads?**

## Tiny Docs 
 The class is initialized by: 
```
import sam_reader
my_files = sam_reader('my_results_folder/')
```
 Def | Input| Return| Kwargs/Other 
--|--|--|--
 sam_to_bam() *static method* | *string* - path to a SAM file | *string* - path to the newly converted and indexed BAM file | Makes system calls to [samtools](http://www.htslib.org/).
 index_bam() *static method* | *string* - path to a BAM file | None | Creates an index file for a BAM file in the same directory using [samtools](http://www.htslib.org/).
 read_counts() *static method* | *string* - path to a BAM file | 
 

## I skim and only look at pictures
![readerhits](/test_files/sam_reader.hits.png)
