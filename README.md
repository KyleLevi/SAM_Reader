
# SAM_Reader

SAM_Reader is a Python package for working with many SAM/BAM formatted files at once. It is built upon [pysam](https://pysam.readthedocs.io/en/latest/api.html) which is great if you are working with a single SAM/BAM file. 

# Quick Start
## Install sam_reader with pip
Open your favorite terminal with pip and type:
```
pip3 install --user sam_reader
```



## Import sam_reader and initialize Sam_Reader:
```
from sam_reader import Sam_Reader
my_files = Sam_Reader('~/pathto/results/')
```
If you are opening SAM files in this directory, they can be automatically converted to BAM files, sorted, and indexed, if you do the following:
```
my_files = Sam_Reader('~/pathto/results/', check_files=True, convert_files=True)
```


## How many reads mapped to each organism in each file?
```
my_files.hits(write_file='new_file.csv')
```
This will loop over every file, and output a CSV file called "new_file.csv" that looks like:

File | Genome  |  Percent Coverage  |  Total Mapped Reads  |  Mapped Reads > 50 bp
--|--|--|--|-- 
SRR3403834 | NC_001416.1 | 99.7 | 6078 | 6024
SRR3403834 | JQ_00887.1 | 33.4 | 5146 | 0
.. | .. | .. | .. | ..


## How do I examine the matches at each individual position in the genome?
```
my_files.per_base_stats(write_file='new_file.csv', organism='my_organism')
```
If you do not specify ```organism=''``` AND there is more than 1 organism present in the files, then you will be shown a list of all organisms and asked to select one at runtime.

Position | Consensus | Percent | A | C | G | T | N | Gap
--|--|--|--|--|--|--|--|--
0 | A | 90.0 | 900 | 83 | 8 | 4 | 5 | 0
1 | C | 100 | 0 | 870 | 0 | 0 | 0 | 0 
.. | .. | .. | .. | ..| .. | .. | .. | ..


## How do I split my files by organism, instead of by dataset?
```
for genome in my_files.genome_lengths.keys():
    my_file.cat('only_' + genome, organism=genome)
```
additionally, you can add ```min_read_len=50``` to only consider matches longer than 50 base pairs.

## I know a little bit of python, is there a way to loop over every read across all files?
Yes! You can use .reads(), which yields one [pysam Aligned Segment object](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment) at a time.  The following example will loop over every read and print the [alignment length](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.query_length) of each.
```
for read in my_files.reads()
    print(read.query_alignment_length())
```
You can also choose to only loop over reads from a specific organism. In this example, only reads that aligned to "NC_001416.1" will have their [query sequence](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.query_sequence) printed.
```
for read in my_files.reads(organism='NC001416.1):
    print(read.query_sequence())
```

## Is there a quick way to...?
Maybe. Open an issue [here](https://github.com/KyleLevi/SAM_Reader/issues). There are probably other people trying to do the same thing, so if it can be done quickly I will absolutely add that code to the docs.


# Tiny Docs 
 The class is initialized by: 
```
from sam_reader import Sam_Reader
my_files = Sam_Reader'my_results_folder/')
```
 Def |  What it does
--|--
reads(**) | Yields one read at a time, over all files. Each read is a [pysam](https://pysam.readthedocs.io/en/latest/index.html) Aligned Segment.  [Here](http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment) is a list of things you can do with an Aligned Segment.
 hits(**)  | Creates a single 2d array (list of lists) from all files with the 5 columns:<br>File - Genome - Percent Coverage - Total Mapped Reads - Mapped Reads > 50 bp
 per_base_stats(**) | Creates a 2d array from the matches overlaying each position in a genome with the columns:<br>Position - Consensus - Percent - A - C - G - T - N - Gap<br>**You will need to specify a single organism if more than one is present in the files.*
  sam_to_bam(*s*)<br>*static method*| Makes system calls to [samtools](http://www.htslib.org/) to convert and index SAM files into the same directory. **This will be performed automatically on any SAM files if they are opened with** ```Sam_Reader('my_folder/', check_files=True, convert_files=True)```
primers(**) | Takes in a single int (your desired primer length) and returns the top 100 most conserved sequences of that length in FASTA format. Conservation of a primer is calculated by *multiplying* the conservation of each individual position in the potential primer.   
  cat(**) | Concatenates reads from all BAM files into a single BAM file. **Kwargs can be used to specify a single organism, or enforce match length requirements**.
  ****kwargs** | Any method with ** can be modified by the following [key word arguments](https://docs.python.org/3/tutorial/controlflow.html#keyword-arguments):<br>organism='my_genome'<br>only_this_file='my_file.bam'<br>min_read_len=50

 
## Here is a diagram of how most defs work.
![readerhits](/test_files/sam_reader.hits.png)

# FAQ
### What do you mean by "Mapped Reads > 50 bp"?
When using .stats(), this column represents the number of matches that aligned more than 50 base pairs. This is to help distinguish strong matches from short matches.
