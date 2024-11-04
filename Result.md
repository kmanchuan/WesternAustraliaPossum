---
share: "true"
---
# Result

## Sequencing and quality control

About 106.29 Gb raw data of sequencing library was generated, containing low sequencing adaptors.

![Adaptor content of R1](../Images/Pasted image 20241104154724.png)

![[Pasted image 20241104155008.png|Adatpor content of R2]]


```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s1.2.Fastqc_out
```

---
After filtering, about 104.12 Gb clean data was obtained.


```bash
#/nesi/project/massey04238/01.possum/00.raw_data/01.WA_wgs_data/s2.filtered
```

---
Insert size of this library was inferred as 400 bp through the intermediate information from Samtools[@danecek2021].  Around 97.37% of clean reads were aligned to reference genome, most of which (~95.92 Gb) were aligned to chromosomes.

| Library     | Raw data (Gb) | Clean data | Mapped rate (%) | Coverage (%) | Mean depth (X) |
| ----------- | ------------- | ---------- | --------------- | ------------ | -------------- |
| 22-M-012721 | 106.29        | 104.12     | 97.37           | 98.94        | 28.82          |

![[Pasted image 20241104165006.png|Statistics of reads mapped to chromosomes]] 

