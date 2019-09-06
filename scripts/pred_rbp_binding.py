import warnings
import os
import time
import shutil
import pandas as pd
import pybedtools
from pybedtools import BedTool, parallel
import kipoi

import argparse

import logging

verbose = False 

def vprint(msg):
    if verbose:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}]\t{msg}")

def predict(bed_file, model, tosample, fasta_file, gtf_file, output_dir):
    vprint(f'========== {bed_file} ==========')
    # bootstrap since too many sequences in single run
    bed = BedTool(bed_file)

    if tosample:
        sample_size = 5000
        iterations = 10
        file_prefix = '.'.join(bed_file.split(".")[:-1]) + f'.{sample_size*iterations}'
        vprint(f'Performing bootstrap sampling (n={sample_size}, iterations={iterations})...')

        # sample
        bed_df = bed.to_dataframe()
        sample_df = bed_df.sample(sample_size*iterations, replace=True)

    # write to file
    intervals_file = f'{file_prefix}.bed'
    
    if tosample:
        BedTool.from_dataframe(sample_df).moveto(intervals_file)
    else:
        intervals_file = bed

    dl_kwargs = dict(
        intervals_file=intervals_file,
        fasta_file=fasta_file,
        gtf_file=gtf_file,
        use_linecache=True,
    )

    vprint('Predicting...')
    pred = model.pipeline.predict(dl_kwargs, batch_size=250)

    out_file = f'{output_dir}/{file_prefix.split("/")[-1]}.{args.model}_pred.txt'
    vprint(f'Saving to {out_file}...')
    
    os.makedirs(output_dir, exist_ok=True)
    with open(out_file, "w") as f:
        f.write('\n'.join([str(_[0]) for _ in pred]))

    vprint('Done.')

if __name__ == '__main__':
    logging.disable(1000)
    warnings.filterwarnings("ignore")
    
    parser = argparse.ArgumentParser(description='Predicting the binding affinity of an RBP to 101bp sequences.')
    parser.add_argument("model", help="RBP model to use to calculate binding affinity")
    parser.add_argument("genome", help="human genome fasta file")
    parser.add_argument("gtf", help="human genome gtf annotation file")
    parser.add_argument("output_dir", help="output directory")
    parser.add_argument("bedfiles", help="input bed file(s) containing 101bp intervals", nargs=argparse.ONE_OR_MORE)
    parser.add_argument("--sample", help="if specified, samples w/ replacement; otherwise attempts to process entire file", default=True, action="store_false")
    parser.add_argument("--output_filename", help="if not specified, defaults to bed filename prefix.samplesize.model_pred.txt (e.g. input=control.bed, output_filename=control.50000.SRSF1_pred.txt)")
    parser.add_argument("--verbose", help="print log messages", action="store_true")
    args = parser.parse_args()

    verbose = args.verbose
    vprint('Verbose printing on.')

    model = kipoi.get_model(f"rbp_eclip/{args.model}")
    
    for bed_file in args.bedfiles:
        predict(bed_file, model, args.sample, args.genome, args.gtf, args.output_dir)
