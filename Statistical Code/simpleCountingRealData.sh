#Find number of regions that intersect: 125266
bedtools intersect -a consindelsHgMmCanFam.bed -b musclehsmm_umbilicalhuvec_stemh1_lungnhlf_skinnhek.bed | wc -l

#Find number of regions in H3K27ac and not in species conservation: 2553039
bedtools subtract -a consindelsHgMmCanFam.bed -b musclehsmm_umbilicalhuvec_stemh1_lungnhlf_skinnhek.bed | wc -l

#Find number of regions in species conservation and not in H3K27ac: 164589
bedtools subtract -a musclehsmm_umbilicalhuvec_stemh1_lungnhlf_skinnhek.bed -b consindelsHgMmCanFam.bed | wc -l