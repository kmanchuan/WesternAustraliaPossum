---
share: "true"
---
# Result

## 1 Sequencing and quality control

### 1.1 Raw data statistics

About 106.29 Gb raw data of sequencing library was generated, containing low sequencing adaptors.

![Adaptor content of R1](./Images/Pasted%20image%2020241104154724.png)

![Adatpor content of R2](./Images/Pasted%20image%2020241104155008.png)

### 1.2 Filtering

After filtering, about 104.12 Gb clean data was obtained.

### 1.3 Mapping to reference genome

Insert size of this library was inferred as 400 bp through the intermediate information from Samtools [@danecek2021].  Around 97.37% of clean reads were aligned to reference genome, most of which (~95.92 Gb) were aligned to chromosomes.

| Library     | Raw data (Gb) | Clean data | Mapped rate (%) | Coverage (%) | Mean depth (X) |
| ----------- | ------------- | ---------- | --------------- | ------------ | -------------- |
| 22-M-012721 | 106.29        | 104.12     | 97.37           | 98.94        | 28.82          |

![Statistics of reads mapped to chromosomes](./Images/Pasted%20image%2020241104165006.png) 

GC-Depth distribution was obtained from mapped reads. GC content concentrated around 40% with solo peak, which reflected the reference genome almost contained no contamination from exotic organisms. A single depth peak was indicated around 30X, which was reasonable according to the data volume of this library. In addition, there was a faint peak near 13X of depth, which might indicate the highly repeat region of the reference genome, such as satellite DNA or large tandem repeats. Because of the limitation of short-read sequencing platform, those complex regions might be difficult to uniquely map, which resulted in low coverage or no coverage in those regions [@liao2019; @pop2008; @alkan2011]. Correspondingly, after checking the sequences near the low-coverage peak manually, plenty of repetitive sequences were observed.

![GC-Depth distribution based on reference genome](./Images/Pasted%20image%2020241105150537.png)

![Repetitive sequences from low-coverage area](./Images/Pasted%20image%2020241106111153.png)

### 1.4 Genome size estimation with *K-mer*

The estimated genome size can be obtained with this fomula: kmer_num/pkdepth, where the kmer_num is the total number of *k-mers*, while pkdepth refers to the most frequent peak of *k-mers*.

Based on kmerfreq result, total number of 17mer was 92,922,327,604 and the main peak of *kmer* frequency was at 27X. The genome size was estimated as 3.44 Gb. Besides the main peak, which could be considered as homo-peak, a repetitive peak was also observed at 54X where the frequency was two times as much as homo-peak. No apparent hetero-peak was obtained, which mean western Australia possums' heterozygosity might be at a low level.

![Estimated genome size by kmerfreq](./Images/Pasted%20image%2020241105163110.png)

By contrast, the estimated genome size based on Jellyfish and Genomescope was between 1.48 Gb with about 5.53% heterozygosity, almost as half as the expected size.
![Estimated genome size by Genomescope](./Images/Pasted%20image%2020241105163838.png)

To confirm whether the result of Genomescope was accurate, genome size was estimated manually based on the result of Jellyfish following this methodology ([Genome Size Estimation Tutorial](https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/)). After filtering low frequency *k-mers* with high numbers (between 1-9) that probably was caused by sequencing error, total 17mer number was 77,623,206,436 and the homo-peak was at 24. As a result, the genome size was estimated to be 3.23 Gb which was close to the estimation by kmerfreq (3.44 Gb) and the actual size of the reference genome (3.36 Gb).

![Estimated genome size based on the result from Jellyfish](./Images/Pasted%20image%2020241106121335.png)

## 2 *De novo* assembly

### 2.1 Assemblers and reduction

For Minia, with *K-mer* size growing, the longer contig N50 and larger genome size was obtained. The assembly with K59 was 3.08 Gb and the contig N50 was 2.79 kb. We didn't test greater *K-mer* size, because compared with Minia under the same *K-mer* size of 59, Megahit could generate longer contig N50 as well as larger total size of genome. So, Megahit was preferred to generate contigs for downstream analysis.

| Assembly       | Total size (Gb) | Contig Number | Contig N50 (kb) | GC (%) |
| -------------- | --------------- | ------------- | --------------- | ------ |
| Minia – K27    | 2.16            | 3,414,553     | 0.81            | 39.47  |
| Minia – K31    | 2.46            | 2,843,843     | 1.23            | 39.40  |
| Minia – K59    | 3.08            | 2,349,099     | 2.79            | 39.39  |
| Megahit – K21  | 0.64            | 2,203,668     | 0.28            | 41.69  |
| Megahit – K29  | 2.32            | 3,243,424     | 1.03            | 39.58  |
| Megahit – K39  | 2.75            | 2,542,296     | 1.78            | 39.69  |
| Megahit – K59  | 3.20            | 2,321,414     | 3.60            | 39.82  |
| Megahit – K69  | 3.32            | 2,124,056     | 3.59            | 39.44  |
| Megahit – K79  | 3.46            | 2,118,584     | 7.06            | 39.58  |
| Megahit – K89  | 3.54            | 1,925,884     | 9.50            | 39.58  |
| Megahit – K99  | 3.40            | 3,875,063     | 1.19            | 39.61  |
| Megahit – K119 | 3.87            | 3,201,053     | 4.79            | 39.52  |
| Megahit – K141 | 4.18            | 3,482,437     | 6.41            | 39.53  |

![Genome size and Contig N50 statistics](./Images/Pasted%20image%2020241106144655.png)

For Megahit, the same trend of contig N50 with Minia was observed between K21 and K89. However, when *K-mer* size exceeded 89, contig N50 dropped rapidly from 9.50 kb in K89 to 1.19 kb in K99, which means the genome sequences were more fragmented and less contiguous. So, we chose the K89 as the best performance.

Given the estimated genome size of western Australia possum was 3.23 -3.44 Gb, the total size of K89 assembly (3.54 Gb) was greater than expectation, which might be caused by redundant sequences. So, reduction process was performed on K89 assembly by Redundans and Purge_dups. After removing redundant sequences, 3.32 Gb genome was obtained by Redundans, while 3.21 Gb was obtained by Purge_dups, and both of the two sizes were close to estimated genome size of western Australia possum.

### 2.2 Scaffolders

Because both Minia and Megahit only generate contigs, to obtain longer and more contiguous genome region, scaffolding process should be performed on those contigs.

To evaluate which scaffolders could performs better with our data,  we firstly compared BESST with Soapscaff based on Megahit - K59 assembly. On the one hand, BESST with default parameters generated the longer scaffold N50, reaching 15.06 kb, and increasing the number of iteration did not seem to improve the result, as scaffold N50 was the same with default parameters version and gap size was even increased by 0.01 Mb. By contrast, the scaffold N50 of the result with *-m 400* parameters was decreased slightly to 14.62 kb, but its gap size dropped dramatically from 85.94 Mb to 43.72 Mb. On the other hand, Soapscaff generated lower scaffold N50 (11.37 kb) as well as smaller gap size (62.03 Mb) compared with BESST with default parameters. In short, Soapscaff could achieve shorter scaffold N50 with relatively less gap, meanwhile BESST could reach the longest scaffold N50 with default parameters, and trade off the length of scaffold N50 for less gap content.
 
| Scaffolder            | Total size (Gb) | Scaf num  | Scaf N50 (kb) | Gap (Mb) | GC (%) |
| --------------------- | --------------- | --------- | ------------- | -------- | ------ |
| BESST (default)       | 3.26            | 1,566,393 | 15.06         | 85.94    | 39.82  |
| BESST (-m 400)        | 3.22            | 1,562,263 | 14.62         | 43.72    | 39.82  |
| BESST (-iter 1000000) | 3.26            | 1,566,352 | 15.06         | 85.95    | 39.82  |
| SoapScaff (-K 41)     | 3.22            | 1,299,394 | 11.37         | 62.03    | 39.81  |

Next, we compared BESST and Redundans based on Megahit - K89 assembly. Both of the two programs were run under default parameters. Unfortunately, scaffolding module of Redundans didn't work well on our data, showing no any improvement after its process.
  
| Scaffolder      | Total size (Gb) | Scaf num  | Scaf N50 (Kb) | Gap (Mb) | GC (%) |
| --------------- | --------------- | --------- | ------------- | -------- | ------ |
| BESST (default) | 3.57            | 1,467,147 | 28.81         | 50.38    | 39.58  |
| Redundans       | 3.54            | 1,925,884 | 9.50          | 0        | 39.58  |

To check whether Megahit contigs also outperformed Minia contigs in scaffolding, we compared their contigs assemblies in K59. As a result, the scaffolds N50 based on Minia contigs was 16.30 kb which was longer 1.24 kb than that of Megahit, but its gap size (198.74 Mb) was larger than twice as much as that of Megahit (85.94 Mb). In addition to that Minia usually generates shorter contig N50 than Megahit, as we described in *Assemblers and reduction* section, we comprehensively chose Megahit as our primary assemblers.

| Assembler     | Scaffolder      | Total size (Gb) | Scaf num  | Scaf N50 (kb) | Gap (Mb) | GC (%) |
| ------------- | --------------- | --------------- | --------- | ------------- | -------- | ------ |
| Minia - K59   | BESST (default) | 3.26            | 1,395,544 | 16.30         | 198.74   | 39.38  |
| Megahit - K59 | BESST (default) | 3.26            | 1,566,393 | 15.06         | 85.94    | 39.82  |

In short, after these comparisons between different assemblers, scaffolders, and settings within the same scaffolders, Megahit might be our optimal assembler, while BESST with *-m 400* parameters might a better option of scaffolder, although slight shorter scaffold N50 often comes together.

Currently, the best version of scaffolds assembly was from BESST with *-m 400*  parameters based on Megahit - K89 contigs after reduction by Purge_dups.

| Scaffolder     | Total size (Gb) | Scaf num | Scaf N50 (Kb) | Gap (Mb) | GC (%) |
| -------------- | --------------- | -------- | ------------- | -------- | ------ |
| BESST (-m 400) | 3.23            | 554,852  | 30.19         | 28.96    | 39.53  |
### 2.3 Gapclosers

For GMcloser, we attempted to use Minia - K59 contigs to fill the gaps of the scaffolds from Megahit and BESST, as different assemblers might have distinct assembling bias. So, contigs from different assemblers may presumably complement each other. However, none of submitted jobs result was completed successfully from GMcloser by now, even though we have adjusted various parameters.

For Redundans, its gapclosing module was tested along during scaffolding. As we described in *Scaffolders* section, its scaffolding module, as well as gapclosing module, didn't work well on our data. So, we discontinue this program for gapclosing.

Gapcloser was robust to close the gaps from both scaffolds generated by its internal pipeline, namely Soapscaff, and external scaffolds like BESST. After gapclosing, the gap size of Soapscaff scaffolds dropped greatly from 62.03 Mb to 7.09 Mb, while contig N50 rose from 3.61 kb to 5.99 kb. Correspondingly, the scaffolds from BESST (based on Megahit - K89 after reduction by Purge_dups) were also improved drastically as their gap size decreased by 73.57% and the contig N50 more than doubled, up to 22.35 kb.

| Version           | Stage     | Total size (Gb) | Gap (Mb) | Scaffold N50 (kb) | Contig N50 (kb) |
| ----------------- | --------- | --------------- | -------- | ----------------- | --------------- |
| Megahit – K59     | contig    | -               | -        | -                 | 3.60            |
| +Soapscaff        | scaffold  | 3.22            | 62.03    | 11.37             | 3.61            |
| ++Gapclolser      | gapfilled | 3.20            | 7.09     | 11.27             | 5.99            |
| Megahit – K89     | contig    | -               | -        | -                 | 9.50            |
| +Purge_dups       | contig    | -               | -        | -                 | 11.04           |
| ++BESST (default) | scaffold  | 3.25            | 51.54    | 31.88             | 11.06           |
| +++Gapclolser     | gapfilled | 3.23            | 13.62    | 31.61             | 22.35           |

By now, the result of gapclosing based on BESST with *-m 400* parameters hasn't obtained. 

So, the optimal *De novo* pipeline for our data was Megahit as an assembler, Purge_dups as a redundant sequences remover, BESST as a scaffolder, and Gapcloser as a gapcloser. Otherwise, given Gapcloser was so potent to reduce the gap size in scaffolds from Soapscaff, decreased by 88.57%, and they are both incorporated in Soapdenovo2, there might be some internal optimizations between them. In other words, an alternative pipeline of Soapscaff-to-Gapcloser maybe generate better result. And we are verifying this thought, waiting for the result.

![Optimal pipeline of *De novo* assembly](./Images/Pasted%20image%2020241107162949.png)

| Version          | Stage     | Total size (Gb) | Gap (Mb) | Scaffold N50 (kb) | Contig N50 (kb) |
| ---------------- | --------- | --------------- | -------- | ----------------- | --------------- |
| Megahit – K89    | Contig    | -               | -        | -                 | 9.50            |
| +Purge_dups      | Contig    | -               | -        | -                 | 11.04           |
| ++BESST (-m 400) | Scaffold  | 3.23            | 28.96    | 30.19             | 11.06           |
| +++Gapclolser    | Gapfilled | Waiting         |          |                   |                 |
| ++Soapscaff      | Scaffold  | Waiting         |          |                   |                 |
| +++Gapcloser     | Gapfilled | Waiting         |          |                   |                 |

### 2.4 GC-Depth distribution

According to the scatter plots, GC-Depth distribution for each assembly was similar with that of reference genome, which indicated that both Minia and Megahit were relatively reliable assemblers, generating contigs without apparent errors. In addition, a low-coverage peak was also observed in each GC-Depth plots. As we explained in *Mapping to reference genome* section, those regions were supposed to be repetitive structure in genome, which were hard to assemble for short reads. And the low-coverage area became deeper along with the increase of *K-mer*, which reflected a fact that relatively greater *K-mer*  could help assembliers solve complex regions in genome [@liao2019].

![GC-Depth distribution of *De novo* assemblies](./Images/Pasted%20image%2020241108102839.png)

### 2.5 Genome completeness

Firstly, BuscoV5 was used to assess the completeness of our *De novo* genomes. However, it seemed to have some conflicts with NeSI platform ([NeSI](https://www.nesi.org.nz/)) where all our jobs were performed, according to NeSI's support team. Following their instruction, we used Compleasm to replace BuscoV5, meanwhile BuscoV4 was also employed to make assessment.

To verify the difference between BuscoV4 and BuscoV5, Megahit - K59 assembly was assessed respectively by each version. The result showed that the difference was very small. The proportion of complete genes of BuscoV5 was 44.8%, slightly lower than that of BsucoV4 at 45.5%. The difference in the proportion of fragmented genes was only 0.7%, and the proportion of missing genes was the same. In conclusion, the difference between the two versions of Busco was negligible.

By contrast, according to the results tested on Minia - K59 scaffolds and Megahit - K89 scaffolds, Compleasm appeared to obtain roughly 6.3 percentage points more complete genes than BuscoV4, as well as approximately 5.7 percentage points more fragmented genes and 12.1percentage points less missing genes. However, Compleasm runs faster than Busco, especially on large and chromosome-level genomes.

For both contigs and scaffolds assessment, the graph showed the larger *K-mer* was set during assembling, the more completeness the assemblies could reach. After reduction, more complete genes were found than before, even though the proportion increased by only 0.1 - 0.2 percent, which means our reduction process was reliable. BESST appeared better completeness than Soapscaff, as the proportion of complete genes in Megahit - K59 scaffolds by BESST was 44.8% - 45.5%, but only 40.9% in scaffolds by Soapscaff. Furthermore, Gapcloser didn't significantly improve completeness, which only increased by 1 percentage point.

![Genome completeness assessment](./Images/Pasted%20image%2020241108131452.png)

### 2.6 Heterozygosity and mutation rate

After calculating heterozygosity or mutation rate of each selected gene, average heterozygosity and mutation rate was obtained. The heterozygosity within western Australia possum was 0.00097, which was at a very low level. On the other hand, the mutation rate between western and eastern Australia possums was 0.0059, which might indicate that the difference between western and eastern Australia possums was larger than the difference in western Australia possum. However, it's noteworthy that the result based on this method was far from being accurate, as we didn't strictly qualify reads alignment or variants calling, not to mention more population data should be involved. This estimation was just a quick look about the complexity of the western Australia possum genome. 

![Average heterozygosity and mutation rate](./Images/Pasted%20image%2020241108144733.png)

## 3 Reference-guided assembly

### 3.1 Consensus

The statistics of the consensus seemed almost the same with reference genome. Especially, the gap size in the consensus and reference genome was identical, 18,696,463 bp, and scaffold number and contig number are also the same. So we suspected that BCFtools might just directly modify the reference genome referring to the VCF file and introduce reference sequence into the consensus in the regions that were not mapped by reads.

| Assembly      | Total size (Gb) | Scaffold Number | Scaffold N50 (Mb) | Gap size (bp) | Contig Number | Contig N50 (Mb) | GC (%) |
| ------------- | --------------- | --------------- | ----------------- | ------------- | ------------- | --------------- | ------ |
| BCF_consensus | 3.36            | 211             | 442.4             | 18,696,463    | 2,007         | 4.31            | 39.48  |
| Public genome | 3.36            | 211             | 442.6             | 18,696,463    | 2,007         | 4.31            | 39.28  |

We manually checked the alignment of NW_023494386.1, and observed that BCFtools did introduce those regions that were not mapped by reads from the reference into the consensus. In that case, the consensus from BCFtools was not credible.

![Checking the consensus](./Images/Pasted%20image%2020241108155530.png)

### 3.2 Subgroups assembly

For now, clean reads are aligning to chromosome 1 of the reference genome. There is no enough result to be discussed yet.

### 3.3 Reference-guided scaffolding

Ragtag scaffolds seemed good as its scaffold N50 reached incredibly 463.37 Mb, while Redundans only achieved 11.55 kb, even lower than that of BESST scaffolds without using reference.

| Scaffolder      | Total size (Gb) | Scaf num  | Scaf N50 (kb) | Gap (Mb) | GC (%) |
| --------------- | --------------- | --------- | ------------- | -------- | ------ |
| Ragtag          | 3.81            | 828,362   | 463.37 (Mb)   | 264.72   | 39.58  |
| Redundans       | 3.92            | 1,901,386 | 11.55         | 377.85   | 39.58  |
| ntJoin          | Waiting         |           |               |          |        |
| BESST (default) | 3.57            | 1,467,147 | 28.81         | 50.38    | 39.58  |

When we checked the AGP file generated by Ragtag, many fixed-length gap were inserted between contigs. Although Ragtag could infer precise length of gap during scaffolding, when that inferring failed, it would insert a fixed-length gap (which could be set by *-g* parameters) between adjacent contigs in the same way with Hi-C assembly. So, the gap size of Ragtag scaffolds wasn't accurate.

![AGP file of Ragtag](./Images/Pasted%20image%2020241108170921.png)