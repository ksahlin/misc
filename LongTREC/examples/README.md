### Run instructions

```
minimap2 --eqx -a [PARAM] genome.fa gene.fa > mm2.sam
python evaluate_instance.py --gtf annotation.gtf  --samfile mm2.sam
```

### Output 

The output consists of four tab-separated values, e.g. `#E	#A	E_min	A_min`, where:

* #E: Number of exact exon alignments.
* #A: Number of exact OR approximate exon alignments (alnmt overlaps true exon).
* E_min: Smallest exon with an exact exon alignment.
* A_min: Smallest exon with an exact OR approximate exon alignment.


### Datasets


* dataset1: no errors, no canonical splices AG-GT
* dataset2: 7% errors, no canonical splices AG-GT
* dataset3: no errors, only canonical splices AG-GT
* dataset4: 7% errors, only canonical splices AG-GT
