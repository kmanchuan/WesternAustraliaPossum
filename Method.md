---
share: "true"
---
# Method
## 1 Sequencing and quality control

### 1.1 Raw data statistics 

Paired-end library with a 400 bp insert size was sequenced on Illumina NovaSeq platform by AGRF, generating 150 bp pared-end reads. Adaptor sequences were identified by FastQC (v0.12.1) [@andrews2010].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s1.2.Fastqc.sh
module load FastQC/0.12.1
fastqc -t 8 -o s1.2.Fastqc_out 22-M-012721_H5L5VDSX7_CTCCACTAAT-AACAAGTACA_L001_R1.fastq.gz 22-M-012721_H5L5VDSX7_CTCCACTAAT-AACAAGTACA_L001_R2.fastq.gz
```

---
### 1.2 Filtering

Raw reads with the rate of *N* higher than 5%, the low-quality bases (quality score < 10) higher than 30%, and the duplicated reads were filtered by SOAPnuke (v 2.1.9) [@chen2018].

```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filter.sh
SOAPnuke filter -1 WGS_raw_R1.fq.gz -2 WGS_raw_R2.fq.gz  -C WGS_clean_R1.fq.gz -D WGS_clean_R2.fq.gz -l 10 -q 0.3 -n 0.05 -c filter_other.config -o s2.filtered -T 10
```

---
### 1.3 Mapping to reference genome

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

Statistics of mapped reads were performed using Bamtools (v2.5.1) [@barnett2011], followed by calculating average coverage of reference genome split by windows of 5,000 bp using BamDeal (v0.27, [GitHub](https://github.com/BGI-shenzhen/BamDeal)). GC-Depth distribution was plotted using in-house scripts, meanwhile specific area sequences near the low-coverage peak (where GC content was between 37.5%-50% while depth was between 10X-17X) were extracted using Seqkit (v2.8.2) [@shen2016].

```shell
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa/s2.stat.sh
bamtools stats -in wa_ref_sorted.bam >wa_ref_sorted.bam.MapRate
bamdeal statistics Coverage -i wa_ref_sorted.bam -r possum_ref.fa -o wa_ref_sorted.bam.stat -w 5000
perl /nesi/project/massey04238/bin/get_wind_tab.pl wa_ref_sorted.bam.stat.DepthGC.gz >win5k.out
Rscript /nesi/project/massey04238/bin/GC_DEPTH.R win5k.out win5k.pdf
```

```shell
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s3.bwa/low.sh
less win5k.out|awk '$5>=37.5&&$5<=50&&$4>=10&&$4<=17'|cut -f1-3|perl -ne 'chomp;my@a=split /\t/,$_;my$s=$a[1]-1;my$e=$a[2];print "$a[0]\t$s\t$e\n"' >low_dep.bed
seqkit subseq --bed low_dep.bed possum_ref.fa -o low_dep.fa
```

---
### 1.4 Genome size estimation with *K-mer*

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

## 2 *De novo* assembly

### 2.1 Assemblers and reduction

Many assemblers were implemented to generate *De novo* assemblies, including Soapdenovo2 (r242) [@luo2012], Abyss (v2.2.5) [@jackman2017], Velvet (v1.2.10) [@zerbino2010], Meraculous (v2.2.6) [@chapmanMeraculousNovoGenome2011], Minia (v0.0.102) [@chikhi2013], and Megahit (v1.2.9) [@li2016]. Among them, only Minia and Megahit worked successfully or generated reasonable result with our data. For example, Abyss and Velvet loaded data for more than 16 hours without any informative response, while Soapdenovo2 just generated an assembly of smaller than 10 Mb, the size of which was obviously incorrect.

To obtain a better assembly, a range of *K-mer* sizes were performed. K27, K31, and K59 were selected in Minia, and K21, K29, K39, K59, K 69, K79, K89, K99, K119, and K141 were chosen in Megahit.

```shell
#/nesi/project/massey04238/01.possum/03.assemble/04.Minia/s1.minia.sh
/nesi/project/massey04238/bin/minia-v0.0.102-bin-Linux/bin/minia -in fq.list -max-memory 100000 -out-tmp /nesi/nobackup/massey04238/temp -nb-cores 10 -out k59_out -kmer-size 59 -out-compress 5
```

```shell
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s1.megahit.sh
/nesi/project/massey04238/bin/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R1.fq.gz -2 /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R2.fq.gz --k-list 99,119,141 --tmp-dir /nesi/nobackup/massey04238/temp -o test_3 -t 10 --continue
```

To remove redundant sequences of genome assemblies, Redundans (v2.01) [@pryszcz2016] and Purge_dups (v1.2.5) [@guan2020] were implemented respectively. The result of Redundans was refined consequently by filtering those potential sequences with both overlap rate and alignment identity more than 90%.

```shell
#Reduans - reduction
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s9.redundans/reduction/redundans.sh
/home/XiaocmFk3eJ/bin/miniconda3/envs/redundans/bin/redundans.py -i /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R*.fq.gz --noscaffolding --nogapclosing -f final.contigs.fa -o redundans -t 10

#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s9.redundans/reduction/redundans/stat.sh
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s9.redundans/reduction/redundans/contigs.reduced.fa.hetero.tsv
less contigs.reduced.fa.hetero.tsv|grep -v "^#"|awk '$4>=0.9&&$5>=0.9'|fish -ff fasta -except - contigs.fa >contigs_reducted_id0.9_cov0.9.fa
```

```shell
#Purge_dups
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/purge_dup.sh
fq1=/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R1.fq.gz
fq2=/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R2.fq.gz

#Index
/opt/nesi/CS400_centos7_bdw/BWA/0.7.18-GCC-12.3.0/bin/bwa index -p K89.contigs final.contigs.fa

#Align
module load SAMtools/1.19-GCC-12.3.0/
bwa mem -t 8 K89.contigs $fq1 $fq2 | samtools view -b -@ 8 -o - > K89.contigs.bam


#purge
/home/XiaocmFk3eJ/bin/purge_dups/bin/ngscstat K89.contigs.bam .
/home/XiaocmFk3eJ/bin/purge_dups/bin/calcuts TX.stat 1>cutoffs 2>calcults.log
/home/XiaocmFk3eJ/bin/purge_dups/bin/split_fa final.contigs.fa >K89.contigs.split
minimap2 -xasm5 -DP -t 8 K89.contigs.split K89.contigs.split | gzip -c - > K89.contigs.split.self.paf.gz
/home/XiaocmFk3eJ/bin/purge_dups/bin/purge_dups -2 -T cutoffs -c TX.base.cov K89.contigs.split.self.paf.gz 1> dups.bed 2> purge_dups.log
/home/XiaocmFk3eJ/bin/purge_dups/bin/get_seqs -e dups.bed final.contigs.fa
```

### 2.2 Scaffolders

Three scaffolding programs were tested, including scaff module of Soapdenovo2 (r242) (Soapscaff), BESST (v2.2.4) [@sahlin2014], and Redundans which actually incorporates SSPACE (v3.0) to perform scaffolding.

```shell
#Soapscaff
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s5.soap-fusion/s1.fushion.sh
SOAPdenovo-fusion -D -K 41 -c final.contigs.fa -g test -p 10 -s final.contigs.config
SOAPdenovo-127mer map -s final.contigs.config -p 10 -g test
SOAPdenovo-127mer scaff -p 10 -g test
```

```shell
#BESST
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s4.besst/k89/run.sh
conda activate besst
/nesi/project/massey04238/bin/BESST/runBESST -c final.contigs.fa -f K89.contigs.sorted.bam -o ./output -orientation fr
```

```shell
#Redundans - scaffolding and gapclosing
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s9.redundans/scaffolding/nonRef_scaff_gapclose/nonRef_scaff_gapclose.sh
conda activate redundans
/home/XiaocmFk3eJ/bin/miniconda3/envs/redundans/bin/redundans.py -i /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R*.fq.gz --noreduction -f final.contigs.fa -o nonRef_scaff_gapclose -t 10 --resume
```

### 2.3 Gapclosers

The gaps of scaffolds can be filled by gapcloser programs using short-reads data. We tested three gapclosers: GMcloser (v1.6.2) [@kosugi2015], Redundans, and Gapcloser (v1.12) [@li2010].

```shell
#GMcloser
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s6.gmcloser/gmcloser/gmcloser.sh
/nesi/project/massey04238/bin/GMcloser-1.6.2/gmcloser -t Scaffolds_pass1.fa  -q Minia_k59_contigs.fa -p test.GMcloser -n 20 -r /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R1.fq.gz /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R2.fq.gz -l 150 -i 400 -d 160 -c
```

```shell
#Redundans - scaffolding and gapclosing
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s9.redundans/scaffolding/nonRef_scaff_gapclose/nonRef_scaff_gapclose.sh
/home/XiaocmFk3eJ/bin/miniconda3/envs/redundans/bin/redundans.py -i /nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered/WGS_clean_R*.fq.gz --noreduction -f final.contigs.fa -o nonRef_scaff_gapclose -t 10 --resume
```

```shell
#Gapcloser
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/s3.besst/k89_purged_besst/m400/gapcloser/gapcloser.sh
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/s3.besst/k89_purged_besst/m400/gapcloser/final.contigs.config
/nesi/project/massey04238/bin/miniconda3/envs/gapcloser/bin/GapCloser -a Scaffolds_pass1.fa -b final.contigs.config -o besstM.gapclosed.fa -l 150 -t 10
```

### 2.4 GC-Depth distribution

Clean reads were mapped and GC-Depth distribution was plotted in the same way as we described in *Mapping to reference genome* section, except that the reference genomes used were our *De novo* assemblies. They were Minia - K31, Megahit - K59, Megahit - K89, and Megahit - K89 after reduction by Purge_dups.

### 2.5 Genome completeness

Busco (v4.1.4, BuscoV4) [@manni2021], Busco (v5.7.1, BuscoV5), and Compleasm (v0.2.6) [@huang2023] were implemented to assess genome completeness using mammalia_odb10 database.

```shell
#BuscoV4
#/scale_wlg_nobackup/filesets/nobackup/massey04238/temp/k89_purged_besstM_busco/buscoV4_sub.sh
/nesi/project/massey04238/bin/miniconda3/envs/buscov4/bin/busco -i Scaffolds_pass1.fa -l /nesi/nobackup/massey04238/temp/busco_database/busco/mammalia_odb10 -o busco_out -m genome -c 20 --offline
```

```shell
#BuscoV5
#/nesi/project/massey04238/01.possum/03.assemble/04.Minia/s2.busco/k31_contig/busco.sh
busco -i test_out.contigs_lenth200.fa -l /home/XiaocmFk3eJ/bin/busco_lineagv5/busco_downloads/lineages/mammalia_odb10/ -o k31.contig.busco -m genome -c 10 --offline --augustus
```

```shell
#Compleasm
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/s4.busco_purged/purged_besstM/k89_purged_besstM_busco/compleasm/compleasm.sh
compleasm run -a Scaffolds_pass1.fa -o compleasm_output -l eukaryota -t 10 -l mammalia -L /nesi/nobackup/massey04238/temp/busco_database/mammalia_odb10 -m busco
```

### 2.6 Heterozygosity and mutation rate

To quick evaluate the heterozygosity of western Australia possum, the first 100 complete and single genes were selected from mammalia_db10. 

```shell
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/s2.heterozygosity/s2.run.sh
grep Complete full_table.tsv |head -100|cut -f3,4,5|perl -ne 'chomp;my@a=split /\t/,$_;my$s=$a[1]-1;print "$a[0]\t$s\t$a[2]\n"' >busco_complete_head100.bed
```

After mapping reads aganist Megahit - K89 contigs after reduction by Purge_dups, variants were called using BCFtools (v1.21) [@danecek2021].

```shell
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/s2.heterozygosity/s1.bcftools_vcf/k89_purged_vcf/run.sh
bcftools mpileup -Ou -f purged.fa wa_purged_sorted.bam | bcftools call -mv -Oz -o output.vcf.gz

bcftools index output.vcf.gz

bcftools norm -f purged.fa output.vcf.gz -Oz -o output.norm.vcf.gz
vcftools --gzvcf output.norm.vcf.gz --recode --recode-INFO-all --stdout --remove-indels|gzip >output.norm.snp.vcf.gz
```

Then, Single Nucleotide Polymorphisms (SNPs) in gene region were extracted using BEDtools (v2.31.0) [@quinlan2010] and the corresponding genotype (GT filed in VCF files) information was split using in-house script. Heterozygosity was calculated roughly by dividing the total SNPs number by the gene length.

```shell
#/nesi/project/massey04238/01.possum/03.assemble/05.megahit/s10.purge_dups/s2.heterozygosity/s3.stat.sh
cat  busco_complete_head100.bed|while read i
        do
        array=(`echo $i|tr ',' ' '`)
        Len=$[${array[2]} - ${array[1]}]
        echo $i|sed 's/ /\t/g' >${array[0]}.bed
        echo "${array[0]}"
        bedtools intersect -a output.norm.snp.vcf.gz -b ${array[0]}.bed > ${array[0]}.snp.out
        less ${array[0]}.snp.out |cut -f10|sed 's/:/\t/g'|cut -f1|sort -k1n|uniq -c |perl -e 'my$out;while(<>){chomp;my@a=split /\s+/,$_;$out.=$a[2]."_".$a[1].";";};print "'${array[0]}'\t'$Len'\t$out\n"' >>busco_complete_head100.snp.stat
        rm ${array[0]}.bed ${array[0]}.snp.out
done
```

In order to estimate the mutation rate between western Australia possums (samples involved in this study) and eastern Australia possums (species of the public genome), short reads of our samples were aligned to the reference, followed by calling variants using BCFtools.

```shell
#/nesi/project/massey04238/01.possum/03.assemble/08.bcftools_consensus/bcftools/bcf_consensus.sh
bcftools mpileup -Ou -f ref.fa wa_ref_sorted.bam | bcftools call -mv -Oz -o output.vcf.gz

bcftools index output.vcf.gz

bcftools norm -f ref.fa output.vcf.gz -Oz -o output.norm.vcf.gz
vcftools --gzvcf output.norm.vcf.gz --recode --recode-INFO-all --stdout --remove-indels >output.norm.snp.vcf
```

Those SNPs within gene region of the first 100 complete and single genes in mammalia_db10 were collected and classified using BEDtools and in-house script. Then, mutation rate was calculated roughly by dividing the total SNPs number by the gene length.

```shell
#/nesi/project/massey04238/01.possum/01.ref/busco/refbusco/s2.heterozygosity/t.sh
cat  busco_complete_head100.bed|while read i
        do
        array=(`echo $i|tr ',' ' '`)
        Len=$[${array[2]} - ${array[1]}]
        echo $i|sed 's/ /\t/g' >${array[3]}.bed
        echo "${array[3]}"
        bedtools intersect -a output.norm.snp.vcf.gz -b ${array[3]}.bed > ${array[3]}.snp.out
        less ${array[3]}.snp.out |cut -f10|sed 's/:/\t/g'|cut -f1|sort -k1n|uniq -c |perl -e 'my$out;while(<>){chomp;my@a=split /\s+/,$_;$out.=$a[2]."_".$a[1].";";};print "'${array[3]}'\t'$Len'\t$out\n"' >>busco_complete_head100.snp.stat
        rm ${array[3]}.bed ${array[3]}.snp.out
done
```

