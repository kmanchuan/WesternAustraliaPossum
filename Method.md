---
share: "true"
---
# Method
## Sequencing and quality control

### Raw data statistics 

Paired-end library with a 400 bp insert size was sequenced on Illumina NovaSeq platform by AGRF, generating 150 bp pared-end reads. Adaptor sequences were identified by FastQC (v0.12.1) [@andrews2010].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s1.2.Fastqc.sh
module load FastQC/0.12.1
fastqc -t 8 -o s1.2.Fastqc_out 22-M-012721_H5L5VDSX7_CTCCACTAAT-AACAAGTACA_L001_R1.fastq.gz 22-M-012721_H5L5VDSX7_CTCCACTAAT-AACAAGTACA_L001_R2.fastq.gz
```

---
### Filtering

Raw reads with the rate of *N* higher than 5%, the low-quality bases (quality score < 10) higher than 30%, and the duplicated reads were filtered by SOAPnuke (v 2.1.9) [@chen2018].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filter.sh
SOAPnuke filter -1 WGS_raw_R1.fq.gz -2 WGS_raw_R2.fq.gz  -C WGS_clean_R1.fq.gz -D WGS_clean_R2.fq.gz -l 10 -q 0.3 -n 0.05 -c filter_other.config -o s2.filtered -T 10
```

---
### Mapping to reference genome

To infer the insert size of the library, the published reference genome from western Australia possums ([GCF_011100635.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011100635.1/)) [@bondAdmixedBrushtailPossum2023] was downloaded from NCBI. Clean reads were aligned against the reference using BWA (v0.7.18-r1243) [@li2013], followed by sorting Bam files using Samtools (v1.19) [@danecek2021].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa/s1.bwa.sh
#Index
/opt/nesi/CS400_centos7_bdw/BWA/0.7.18-GCC-12.3.0/bin/bwa index -p possum_ref possum_ref.fa
#Align
module load SAMtools/1.19-GCC-12.3.0/
fq1=/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R1.fq.gz
fq2=/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R2.fq.gz
ref=possum_ref
/opt/nesi/CS400_centos7_bdw/BWA/0.7.18-GCC-12.3.0/bin/bwa mem -t 20 $ref $fq1 $fq2|samtools view -bS -@ 20 -|samtools sort -@ 20 - -o wa_ref_sorted.bam
```

Statistics of mapped reads were performed using Bamtools (v2.5.1) [@barnett2011], followed by calculating average coverage of reference genome split by windows of 5,000 bp using BamDeal (v0.27, [GitHub](https://github.com/BGI-shenzhen/BamDeal)). GC-Depth distribution was plotted using in-house scripts.


```shell
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa
bamtools stats -in wa_ref_sorted.bam >wa_ref_sorted.bam.MapRate
bamdeal statistics Coverage -i wa_ref_sorted.bam -r possum_ref.fa -o wa_ref_sorted.bam.stat -w 5000
perl /nesi/project/massey04238/bin/get_wind_tab.pl wa_ref_sorted.bam.stat.DepthGC.gz >win5k.out
Rscript /nesi/project/massey04238/bin/GC_DEPTH.R win5k.out win5k.pdf
```

---
### Genome size estimation with *K-mer*

Genome size of western Australia possum was estimated through *K-mer* approach, employing kmerfreq (v1.0) [@liu2013] with 17mer.

```shell
#/nesi/project/massey04238/01.possum/02.jellyfish/kmerfreq/kmer.sh
kmerfreq -k 17 -l fq.list -t 10  >17mer.freq 2>17mer.log
```

Meanwhile, Jellyfish (v2.3.0) [@marcaisFastLockfreeApproach2011] and Genomescope (v1.0) [@vurture2017] were also used to estimate genome size and its heterozygosity.

```shell
#/nesi/project/massey04238/01.possum/02.jellyfish/jellyfish/jellyfish.sh

fq1=/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R1.fq.gz
fq2=/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R2.fq.gz

module load Jellyfish/2.3.0-gimkl-2020a

jellyfish count -C -m 17 -t 20 -s 1000000000 -o 17mer_out <(zcat $fq1) <(zcat $fq2)

jellyfish histo -t 20 17mer_out > 17.histo
Rscript genomescope.R 17.histo 17 150 out 10000
```

## *De novo* assembly

### Assemblers

Many assemblers were implemented to generate *De novo* assemblies, including Soapdenovo2 (r242) [@luo2012], Abyss (v2.2.5) [@jackman2017], Velvet (v1.2.10) [@zerbino2010], Meraculous (v2.2.6) [@chapmanMeraculousNovoGenome2011], Minia (v0.0.102) [@chikhi2013], and Megahit (v1.2.9) [@li2016].