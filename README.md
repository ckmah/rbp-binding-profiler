# rbp-binding-profiler

## Setup
Requires `Python 3.6 and Anaconda`.

1. Install dependencies with file with `pip install -r requirements.txt`
2. Create and activate virtual environment with `kipoi ` env utility.
    - env name **w/out** gpu flag: `kipoi-rbp_eclip`
    - env name **with** gpu flag:  `kipoi-gpu-rbp_eclip`
```
kipoi env create --gpu rbp_eclip source=kipoi
conda activate kipoi-gpu-rbp_eclip
```

## Scripts

- `pred_rbp_binding.py` - main script; predict binding scores for each sequence using "rbp-eCLIP" model from `kipoi` model repository
- `bam2bed.py` - converts ENCODE RBP Knockdown RNA-seq bam files to bed files
