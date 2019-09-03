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
logging.disable(1000)

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("model", help="RBP model to use to calculate binding affinity")
parser.add_argument("bedfiles", help="bed file(s) input", nargs=argparse.ONE_OR_MORE)
parser.add_argument("--verbose", help="print log messages", action="store_true")

args = parser.parse_args()

if args.verbose:
    vprint('Verbose printing on.')


def vprint(msg):
    if args.verbose:
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}]\t{msg}")


project_path = "/home/ubuntu/binding-profiles"
fasta_file = os.path.join(project_path, "data/hg38/hg38.fa")
gtf_file = os.path.join(project_path, "data/hg38/Homo_sapiens.GRCh38.97.gtf")


def predict(bed_file, model):
    vprint(f'========== {bed_file} ==========')
    # bootstrap since too many sequences in single run
    bed = BedTool(bed_file)

    sample_size = 5000
    iterations = 10
    file_prefix = '.'.join(bed_file.split(".")[:2]) + f'.{sample_size*iterations}'
    vprint(f'Performing bootstrap sampling (n={sample_size}, iterations={iterations})...')

    # sample
    bed_df = bed.to_dataframe()
    bootstrap_df = bed_df.sample(sample_size*iterations, replace=True)

    # write to file
    bootstrap_file = f'{file_prefix}.bed'
    BedTool.from_dataframe(bootstrap_df).moveto(bootstrap_file)

    dl_kwargs = dict(
        intervals_file=bootstrap_file,
        fasta_file=fasta_file,
        gtf_file=gtf_file,
        use_linecache=True,
    )

    vprint('Predicting...')
    pred = model.pipeline.predict(dl_kwargs, batch_size=250)

    out_file = f'{project_path}/results/{file_prefix.split("/")[-1]}.{args.model}_pred.txt'
    vprint(f'Saving to {out_file}...')

    with open(out_file, "w") as f:
        f.write('\n'.join([str(_[0]) for _ in pred]))

    vprint('Done.')

model = kipoi.get_model(f"rbp_eclip/{args.model}")

for bed_file in args.bedfiles:
    predict(bed_file, model)
