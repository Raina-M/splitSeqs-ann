#!/bin/python
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-s','--seqsize',help='Length of the sequence going to be split.')
parser.add_argument('-t','--threshold',help='Size limitation for sequences.')
parser.add_argument('-f','--file',help="Annotation file, either in GTF or GFF format.")
parser.add_argument('-o','--outdir',help="Output directory. String.")

args = parser.parse_args()


# find the max gap in every split range
def search_max_gap(df_ann):
	df_pos=df_ann.iloc[:,3:5]
	df_pos.columns=['start','end']
	
	# initialize start and end position
	s=df_pos['start'].iloc[0] 
	e=df_pos['end'].iloc[0]
	gap=0

	for index,row in df_pos.iterrows():
		s_c=row['start'] # current start position
		e_c=row['end'] # current end position
		if s_c>=e:
			gap_c=s_c-e
			if gap_c > gap:
				gap=gap_c
				gap_range=[e, s_c]
				
			s=s_c
			e=e_c
			
		elif e_c>e:
			e=e_c
	
	return gap_range

def compute_split_position(seqsize, sizeLimit, f_ann):

	df_ann = pd.read_csv(f_ann, sep="\t", header=None)
	N=seqsize//sizeLimit
	
	lines=[]
	split_position=0
	for i in np.arange(1,N+1):
		# compute the range of i-th split position
		lower_boundary = seqsize-(N-i+1)*sizeLimit
		upper_boundary = sizeLimit+split_position
		
		# extract the annotation records in the above range to search the max gap
		df_ann_i = df_ann[~((df_ann.iloc[:,4]<lower_boundary) | (df_ann.iloc[:,3]>upper_boundary))]
		
		# find the max gap in i-th split position range
		max_gap_range = search_max_gap(df_ann_i)
		
		# i-th split position is at the midpoint of max gap
		split_position = max_gap_range[0] + (max_gap_range[1]-max_gap_range[0])//2
		
		seqname=df_ann_i.iloc[0,0]
		lines.append([seqname+'_'+str(i), split_position])
	return pd.DataFrame(lines)

def main(args):
	seqsize=int(args.seqsize)
	sizeLimit=int(args.threshold)
	N=seqsize//sizeLimit
	
	# dataframe: 1st col: seqname_N, 2nd col: end position of this sequence segment
	df_split_pos = compute_split_position(seqsize, sizeLimit, args.file)
	
	# start position of every segment
	start_pos = [1]+(df_split_pos.iloc[0:N-1,1]+1).tolist()
	# information of last segment:
	seqname=''.join(df_split_pos.iloc[0,0].split('_')[:-1])
	last_seg = [seqname+'_'+str(N+1), df_split_pos.iloc[N-1,1]+1, seqsize]
	
	# add one column - start position of every segment as the 2nd column
	# then the dataframe has 3 columns:
	# 1st col: seqname_N, 2nd col: start position, 3rd col: end position
	df_split_pos.insert(loc=1, column=None, value=start_pos)
	
	# add last segment to the last line
	df_seg_range=df_split_pos.append(pd.Series(last_seg, index=df_split_pos.columns), ignore_index=True)
	
	# ouput the segment names, start and end postions to a csv format file
	df_seg_range.to_csv(args.outdir, sep="\t", header=None, index=None)
    

if __name__ == "__main__":
    main(args)

