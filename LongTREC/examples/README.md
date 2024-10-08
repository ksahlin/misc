### Requirements

* Have Python installed and able to execute `python` from a terminal window. 
* Have minimap2 installed [minimap2](https://github.com/lh3/minimap2/tree/master?tab=readme-ov-file#getting-started). 


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


* dataset1: no errors, no canonical splices GT-AG
* dataset2: 7% errors, no canonical splices GT-AG
* dataset3: no errors, only canonical splices GT-AG
* dataset4: 7% errors, only canonical splices GT-AG
