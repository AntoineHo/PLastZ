# PLastZ
Easy to use Parallel Lastz script

## Requirements:
- python 3.x (3.7.7 tested)
- Biopython
- samtools
- lastz

## Installation (Bioconda recommended)
```bash
conda create -n plastz
conda activate plastz
conda install biopython lastz samtools
```
Using `requirements.txt`
```bash
conda create -n plastz --file requirements.txt
```

## Quick start
Any alignment:
```bash
python PlastZ.py <query> <target> <outdir> --lastz-options <"quoted string"> --processes <INT>
```

Self-alignments:
```bash
python PlastZ.py <query> <query> <outdir> --lastz-options <"quoted string"> --processes <INT>
```

Exemple:
```bash
python PlastZ.py reads.fa reference.fa read_aligned/plastz_out -p 10 \
  -lo="--format=general:ngap,nmismatch,name1,size1,length1,start1,end1,name2,size2,length2,start2,end2,identity"
```


## Warning:
Depending on the output format, multiple headers may be found in the output single alignement file!

Self alignments with Plastz <query> <query> may result differently from the command `lastz query\[multiple\] --self`

This script only launches many parallel Lastz alignments, using it with many sequences will help improve runtime but is still not as efficient as using another multithreaded software. 
