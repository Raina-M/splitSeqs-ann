# splitSeqs-ann
This script is to split the contigs/sequences whose size is over the pre-defined threshold. The split location is based on the provided corresponding annotation file, specifically, if the size of a certain contig (suppose its name as *contig1*) is over the threshold, it will be split at the midpoint of the largest possible gap between its annotated features, and simultaneously to assure the lengths of all sequences after splitting are all under size limitation. Also, rename the contig as *contig1_1*, *contig1_2*, ..., *contig1_M* (*M* is the number of sequences split from *contig1*) both in sequence file (FASTA) and annotation file (GTF/GFF). The start and end positions in the annotations are properly adjusted based on the split position as well, i.e, subtracting the split site from the position of annotaion features.


Example of application: Running 10X Genomics with a reference contains a contig longer than 536.8Mb (2^29 bp) will not be supported due to limitations of the `.bai` format. In this case, the large contig should be split where it avoids the features. The genome reference fasta and the GTF file will need to be edited together, so that the locations of the features still add up.


```
sh ./split_seqs_and_ann.sh -h
This script is to split the contigs whose size is over a pre-setted value.
The split location is based on the provided corresponding annotation file,
i.e., the contig is splitted at the largest gap between the annotated features.

Syntax: sh ./split_seqs_and_ann.sh -f <FASTA/FA> -a <GTF/GFF> [-s|o|k|h]
options:
-f		  Sequence file that is going to be split.
-a		  Annotation file corresponding to the sequence file.
-s INT		  Length threshold to split the sequences. Default: 2^29nt or 536,870,912nt.
-o DIR		  Output directory. All generated files are put in this directory. 
                  Default: current working directory.
-k		  Keep all intermediate files when specified.
-h		  Print this help and exit.
```

## Dependencies
samtools (not required if index of sequences, i.e., `.fai` file, exists)

python 2.7+


## Input
Mandatory arguments: sequence file (FASTA/FA), annotaions (GTF/GFF)

Options: size limitation, output directory, keep intermediate file or not


## Processing details
1. Find contigs whose size exceeds size limitation, index FASTA -> FAI;
2. Nail down the segmentation positions,
 for every contig:
 - get the quotient of (contig length / size limitation): *Q*, then the contig will be split into *Q+1* segments;
 - compute the segmentation position *P<sub>i</sub>* : len(contig) - (*Q*-*i*+1) x size_limit <= *P<sub>i</sub>* <= size_limit + *P<sub>i</sub>*  - 1, 1 <= *i* < =*Q* and *i* is an integer;
 - find the max gap in every segmentation range by traversal of annotation records in the above range;
 - the split position is at the midpoint of the max gap range.
3. Split annotations, rename contigs, and update start and end positions of split contigs;
4. Split sequence file and rename contigs.


## Output
modified_sequences.fasta, modified_annotaions.gtf/gff

Intermediate files will remain in the output directory when `-k` is activated.
