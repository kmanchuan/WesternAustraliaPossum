---
share: "true"
---
# Method
## Sequencing and quality control

Paired-end library with a 400 bp insert size was sequenced on Illumina NovaSeq platform by AGRF, generating 150 bp pared-end reads. Adaptor sequences were identified by FastQC (v0.12.1) [@andrews2010].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s1.2.Fastqc.sh
module load FastQC/0.12.1
fastqc -t 8 -o s1.2.Fastqc_out 22-M-012721_H5L5VDSX7_CTCCACTAAT-AACAAGTACA_L001_R1.fastq.gz 22-M-012721_H5L5VDSX7_CTCCACTAAT-AACAAGTACA_L001_R2.fastq.gz
```

---
Raw reads with the rate of *N* higher than 5%, the low-quality bases (quality score < 10) higher than 30%, and the duplicated reads were filtered by SOAPnuke (v 2.1.9) [@chen2018].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filter.sh
SOAPnuke filter -1 WGS_raw_R1.fq.gz -2 WGS_raw_R2.fq.gz  -C WGS_clean_R1.fq.gz -D WGS_clean_R2.fq.gz -l 10 -q 0.3 -n 0.05 -c filter_other.config -o s2.filtered -T 10
```

---
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

Statistics of mapped reads were observed using Bamtools (v2.5.1) [@barnett2011], followed by 