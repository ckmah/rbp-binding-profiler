#! /usr/bin/python

import pandas as pd
import swifter
import time
import os
import sys

import pybedtools
from pybedtools import BedTool

from tqdm import tqdm

import requests, json

import warnings
warnings.filterwarnings("ignore")

verbose = True

def vprint(msg):
		if verbose:
				timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
				print(f"[{timestamp}]\t{msg}")


def preprocess_bam_to_bed(bam, output):
    '''
    Given local bam file, convert reads to set of 101bp intervals and output as bed file. Filter for reads thats are 
    '''
		# convert bam to bed
		vprint("Converting bam to bed...")
		bam = BedTool(bam)
		bed = bam.bam_to_bed()

		# filter intervals
		vprint("Filter reads by size...")
		bed_chunk_iter = bed.to_dataframe(chunksize=10000000)  # chunk large file
		chunks = []
		for chunk in bed_chunk_iter:
				keep = (
						chunk[["start", "end"]]
						.swifter.progress_bar(enable=True, desc=bam)
						.apply(lambda row: is_valid_interval(row["start"], row["end"]), axis=1)
				)

				chunks.append(chunk[keep])

		bed_df = pd.concat(chunks)

		# 101bp interval for input
		vprint("Define 101bp intervals...")
		bed_df["end"] = (
				bed_df["start"].swifter.progress_bar(
						enable=True).apply(define_interval)
		)
		bed_df["name"] = "-"

		# remove duplicates
		vprint("Drop duplicate intervals...")
		bed_df.drop_duplicates(inplace=True)

		# TODO extraneous chromosomes?
		vprint("Remove extra chromosomes...")
		chromosomes = list(range(1, 23))
		chromosomes.append('X')
		chromosomes.append('Y')
		chromosomes = [f'chr{c}' for c in chromosomes]
		bed_df = bed_df.loc[bed_df['chrom'].isin(chromosomes)]

		# Save result
		vprint(f"Saving {bed_df.shape[0]} intervals...")
		BedTool.from_dataframe(bed_df).moveto(output)

		# cleanup tmp files
		pybedtools.cleanup(remove_all=True)

		vprint("Done.")


def is_valid_interval(start, end):
    '''Check if read interval is between 90 and 150 bp.'''
		dist = end - start
		return (dist > 90) & (dist < 150)


def define_interval(start):
		return start + 101


def preprocess_wrap(row):
	bam, sample_name = row
	sample_json_url = f'https://www.encodeproject.org/files/{bam}/?format=json'
	sample_json = requests.get(sample_json_url, headers={'accept': 'application/json'}).json()
		
	vprint(f'{bam} - {sample_name}')
	preprocess_bam_to_bed(f'data/bam/{bam}.bam', f'data/bed/{sample_name}.filt.bed')


def main():
		sample_metadata = pd.read_csv(
				'data/bam/rbp_model_metadata.tsv', sep='\t', header=None)
		sample_metadata = sample_metadata.iloc[:,[0,18]]

		sample_metadata.columns = ['bam', 'sample_name']
		sample_metadata['sample_name'] = sample_metadata['sample_name'].apply(lambda x: x.split('-')[0])

		for row in tqdm(sample_metadata.itertuples(name=False, index=False), total=len(sample_metadata), file=sys.stdout):
			preprocess_wrap(row)

if __name__ == '__main__':
		main()
