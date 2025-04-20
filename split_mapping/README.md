# Description

Workflow for analyzing split mappings in a simple simulated genome without repeats.

# Requirements

Requires strobealign, minimap2, and bwa mem in path.

# Running 

Run as follows

```
python simulate_instance.py dataset/ # simulates genome and reads with controlled breakpoints
./run_all.sh dataset/  # runs mm2, bwa and strobealign on sim instance
```

# Pipeline details

### simulate_instance.py

Simulates a contig pair (c1, c2). Each contig is of size 1000nt 
and has a site at B (B=500, 1-indexed for gtf), where a read is split mapped over. 
The first m nucleotides of the read maps to c1, 
and the last R-m nucleotides map on c2. We pick all m in interval [10, R-10], 
resulting in R-10 reads over the same region with a unique splicing. 

The above described experiment is repeated N times (default = 100),
to account for randomness in syncmer/minimizer generation around the region.

While it is a simple experiment, it should serve as a first test dataset to
investigate splice mapping performance. 

### evaluate_instance.py

Derives statistics on the fraction of alingments were correclty split over 
breakpoint with left flank correctly aligned and right flank correctly aligned. 
There are three classes: 
  - exact (exact match to bp)
  - approximate (<5nt from true bp)
  - error (something is aligned but >=5nt from bp)

If the fractions of exact, approximate, and error do not sum to one,
the remaining fraction are missing flanks (nothing is aligned on that side of the bp). 