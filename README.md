# PLastZ
Parallel Lastz

## Installation & Quick Start
```bash
conda create -n plastz
conda activate plastz
conda install biopython lastz samtools
python PlastZ.py <query> <target> <outdir> --lastz-options <"quoted string"> --processes <INT>
```

Self-alignments:
```bash
python PlastZ.py <query> <query> <outdir> --lastz-options <"quoted string"> --processes <INT>
```

Exemple:
```bash
python PlastZ.py reads.fa reference.fa read_aligned/plastz_out -p 10 \
  -lo "--format=general:ngap,nmismatch,name1,size1,length1,start1,end1,name2,size2,length2,start2,end2,identity"
```
## Requirements:
to do...

## Warning:
Results may differ from the command `lastz query\[multiple\] --self`
