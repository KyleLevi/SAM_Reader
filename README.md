
### This repo is still under construction. The code is working, but not totally bug free.

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
 Def |  What it does
--|--
reads(**) | Yields one read at a time, over all files. Each read is a [pysam](https://pysam.readthedocs.io/en/latest/index.html) Aligned Segment.  [Here](http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment) is a list of things you can do with an Aligned Segment.
 hits(**)  | Creates a single 2d array (list of lists) from all files with the 5 columns:<br>File - Genome - Percent Coverage - Total Mapped Reads - Mapped Reads > 50 bp
 per_base_stats(**)| Creates a 2d array from the matches overlaying each position in a genome with the columns:<br>Position - Consensus - Percent - A - C - G - T - N - Gap<br>**You will need to specify a single organism if more than one is present in the files.*
  sam_to_bam(*s*)<br>*static method*| Makes system calls to [samtools](http://www.htslib.org/) to convert and index SAM files into the sam directory. **This will be performed automatically on any SAM files if they are opened with** ```sam_reader('my_folder/', check_files=True, convert_files=True)```
  cat(**) | Concatenates reads from all BAM files into a single BAM file. **Kwargs can be used to specify a single organism, or enforce match length requirements**.
  ** | Any method with ** can be modified by the following [key word arguments](https://docs.python.org/3/tutorial/controlflow.html#keyword-arguments):<br>organism='my_genome'<br>only_this_file='my_file.bam'<br>min_read_len=50

 
## Here is a diagram of how most defs work.
![readerhits](/test_files/sam_reader.hits.png)
