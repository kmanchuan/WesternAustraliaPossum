---
share: "true"
---
# Result

## Sequencing and quality control

###  Raw data statistics

About 106.29 Gb raw data of sequencing library was generated, containing low sequencing adaptors.

![Adaptor content of R1](./Images/Pasted%20image%2020241104154724.png)

![Adatpor content of R2](./Images/Pasted%20image%2020241104155008.png)

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s1.2.Fastqc_out
```

---
### Filtering

After filtering, about 104.12 Gb clean data was obtained.

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered
```

---
### Mapping to reference genome

Insert size of this library was inferred as 400 bp through the intermediate information from Samtools[@danecek2021].  Around 97.37% of clean reads were aligned to reference genome, most of which (~95.92 Gb) were aligned to chromosomes.

| Library     | Raw data (Gb) | Clean data | Mapped rate (%) | Coverage (%) | Mean depth (X) |
| ----------- | ------------- | ---------- | --------------- | ------------ | -------------- |
| 22-M-012721 | 106.29        | 104.12     | 97.37           | 98.94        | 28.82          |

![Statistics of reads mapped to chromosomes](./Images/Pasted%20image%2020241104165006.png) 


GC-Depth distribution was obtained from mapped reads. GC content concentrated around 40% with solo peak, which reflected the reference genome almost contained no contamination from exotic organisms. A single depth peak was observed around 30X, which was reasonable according to the data volume of this library. In addition, there was a faint peak near 13X of depth, presenting the reference genome might contains some heterozygous sequences.

![GC-Depth distribution based on reference genome](./Images/Pasted%20image%2020241105150537.png)


```shell
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa/wa_ref_sorted.bam.MapRate

#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa/win5k.out

#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa/win5k.pdf
```

---
### Genome size estimation with *K-mer*

The estimated genome size can be obtained with this fomula: kmer_num/pkdepth, where the kmer_num is the total number of *k-mers*, while pkdepth refers to the most frequent peak of *k-mers*.

Based on kmerfreq result, total number of 17mer was 92,922,327,604 and the main peak of *kmer* frequency was at 27X. The genome size was estimated as 3.44 Gb. Besides the main peak, which could be considered as homo-peak, a repetitive peak was also observed at 54X where the frequency was two times as much as homo-peak. No apparent hetero-peak was obtained, which mean western Australia possums' heterozygosity might be at a low level.

![Estimated genome size by kmerfreq](./Images/Pasted%20image%2020241105163110.png)

```shell
#/nesi/project/massey04238/01.possum/02.jellyfish/kmerfreq/17mer.freq

#/nesi/project/massey04238/01.possum/02.jellyfish/kmerfreq/17mer.log
```

By contrast, the estimated genome size based on Jellyfish and Genomescope was between 1.48 Gb with about 5.51% heterozygosity.

![Estimated genome size by Genomescope](./Images/Pasted%20image%2020241105163838.png)

```shell
#/nesi/project/massey04238/01.possum/02.jellyfish/jellyfish/out/17.histo

#/nesi/project/massey04238/01.possum/02.jellyfish/jellyfish/out/summary.txt
```

