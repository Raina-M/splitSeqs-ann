#!/bin/bash

# default size limitation
sizeLimit=536870912
# default suffix of split files
suffixID='fitted_size'
# default output path
outDir=$(pwd)
# get the script path
script_full_path=$(dirname "$0")

################################################################################                                                                  
printHelp()
{
   # Display Help
   echo "This script is to split the contigs whose size is over a pre-setted value."
   echo "The split location is based on the provided corresponding annotation file,"
   echo "i.e., the contig is splitted at the largest gap between the annotated features."
   echo
   echo "Syntax: sh ./split_seqs_and_ann.sh -f <FASTA/FA> -a <GTF/GFF> [-s|o|k|h]"
   echo "options:"
   echo "-f		Sequence file that is going to be split."
   echo "-a		Annotation file corresponding to the sequence file."
   echo "-s INT		Length threshold to split the sequences. Default: 2^29nt or 536,870,912nt."
   echo "-o DIR		Output directory. All generated files are put in this directory."
   echo "		Default: current working directory."
   echo "-k		Keep all intermediate files when specified."
   echo "-h		Print this help and exit."
   echo
}

########### HANDLE OPTIONS ##############
while getopts "hf:a:s:o:k" opt;
do
    case $opt in
    f) _fasta=$OPTARG ;;
	a) _ann=$OPTARG ;;
	s) re='^[0-9]+$'
		if [[ $OPTARG =~ $re ]]; then
			echo " -s|--size-limit flag was triggered, will split sequences with size over $OPTARG nt." >&1
			sizeLimit=$OPTARG
		else
			echo ":( Wrong syntax for size threshold parameter value. Using the default value ${sizeLimit}." >&2
			printHelp
			exit 1
		fi ;;
	h) printHelp
		exit ;;
	o) outDir=$OPTARG ;;
	k) echo " --keep-files was triggered, will not remove intermediate files generated." >&1
		keep_file=true ;;
    esac
done


# --------------------- Step 1: Find sequences with size over threshold --------------------
# index sequence file
if [[ ! -f $_fasta.fai ]]
then
    echo "Index of sequence fasta file does not exist on your filesystem, will be created."
    samtools faidx $_fasta
fi

# extract sequence names whose length is over size limitation (default: 2^29nt)
awk -F $'\t' -v sl=$sizeLimit '$2>sl {print $1}' $_fasta.fai > ${outDir}/large_seqs.list

sed '/^#/d' $_ann > ${outDir}/modified_annotations.gtf
touch ${outDir}/large_sequences.fasta

while read -r seqname
do
	[ -d ${outDir}/$seqname ] || mkdir -p ${outDir}/$seqname
	sed '/^#/d' $_ann | awk -F $'\t' -v sn=$seqname '$1==sn {print $0}' > ${outDir}/$seqname/annotation_${seqname}.gtf
	
	# --------------------- Step2: Determine segmentation positions for each large sequences ---------------------
	# for every sequence/contig:
	# - get the quotient of (contig length / size limitation) -> Q, then the contig will be split into Q+1 segments
	# - segmentation position Pi range: len(contig)-(Q-i+1)*size_limit <= Pi <= i*size_limit, 1<=i<=Q and i is an integer
	# - find the max gap in every segmentation range by annotation record traversal in the range
	# - the split position is at the midpoint of the max gap range

	seqsize=$(awk -F $'\t' -v sn=$seqname '$1==sn {print $2}' $_fasta.fai)
	python ${script_full_path}/compute_split_positions.py -s $seqsize  -t $sizeLimit \
							   -f ${outDir}/$seqname/annotation_${seqname}.gtf \
							   -o ${outDir}/$seqname/split_positions.txt
	
	# --------------------- Step 3: Modify annotation file accordingly ---------------------
	# for each sequence/contig with size over threshold:
	# sequence name is changed into original_seqname_N
	# offsets of start & end positions are modified based on the split point/new start position
	
	# NL: number of segments this sequence is split into
	NL=$(wc -l ${outDir}/$seqname/split_positions.txt | cut -d' ' -f1)
	for ((i=1; i<=$NL; i++))
	do
		START=$(awk -v line=$i 'NR==line {print $2}' $outDir/$seqname/split_positions.txt)
		END=$(awk -v line=$i 'NR==line {print $3}' $outDir/$seqname/split_positions.txt)
		echo "$(awk -F'\t' -v sn=$seqname -v snNew=${seqname}_$i -v start=$START -v end=$END \
		        'BEGIN{OFS=FS} { if($1==sn && $4>=start && $5<=end) {$1=snNew; $4=$4-start+1; $5=$5-start+1}} 1' \
		        ${outDir}/modified_annotations.gtf)" > ${outDir}/modified_annotations.gtf
	done

	# --------------------- Step 4: Split sequence files accordingly ---------------------
	# for each sequence/contig with size over threshold:
	
	# 1). extract the sequence from fasta file and remove line breaks
	line_bondary=$(grep -n ">" $_fasta | grep -w -A1 $seqname | cut -d":" -f1)
	start=$(echo $line_bondary | cut -d' ' -f1)
	start_line=$((start+1))
	end=$(echo $line_bondary | cut -d' ' -f2)
	end_line=$((end-1))	
	awk -v s=$start_line -v e=$end_line 'NR >= s && NR <= e' $_fasta \
	| tr -d '\n' > ${outDir}/$seqname/${seqname}.fasta
	 
	# 2). extract sequences based on the start and end positions of each segment
	line_size=$(($(sed -n '2p' $_fasta | wc -c)-1))
	echo "Every line has $line_size nucleotides."
	for ((i=1; i<=$NL; i++))
	do
		START=$(awk -v line=$i 'NR==line {print $2}' ${outDir}/$seqname/split_positions.txt)
		END=$(awk -v line=$i 'NR==line {print $3}' ${outDir}/$seqname/split_positions.txt)
		echo "The start position of segment $i of sequence $seqname is $START"
		echo "The end position of segment $i of sequence $seqname is $END"
		
		# select bases from start to end positions of i-th segment
		# add line break
		# add tag: '>seqname_i'
		cut -c$START-$END ${outDir}/$seqname/${seqname}.fasta | awk '{gsub(/.{'"$line_size"'}/,"&\n")}1' > ${outDir}/$seqname/${seqname}_$i.fasta
		sed -i "1s/^/>${seqname}_${i}\n/" ${outDir}/$seqname/${seqname}_$i.fasta
	done
	
	# 3). merge all segments into a single FASTA file
	cat ${outDir}/$seqname/${seqname}_*.fasta >> ${outDir}/large_sequences.fasta
	# remove intermediate files unless -k flag is activated
	rm ${outDir}/$seqname/${seqname}.fasta
	if [ -z "$keep_file" ]
	then
		rm -r ${outDir}/$seqname
	fi
		
done < ${outDir}/large_seqs.list

# merge all modified and unmodified sequences into a complete FASTA file
awk -v list=${outDir}/large_seqs.list 'BEGIN{while((getline<list)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $_fasta > ${outDir}/normal_sequences.fasta
cat ${outDir}/large_sequences.fasta ${outDir}/normal_sequences.fasta > ${outDir}/modified_sequences.fasta

# remove all intermidiate files unless -k flag is activated
if [ -z "$keep_file" ]
then
    rm ${outDir}/large_seqs.list
    rm ${outDir}/large_sequences.fasta
    rm ${outDir}/normal_sequences.fasta
fi
