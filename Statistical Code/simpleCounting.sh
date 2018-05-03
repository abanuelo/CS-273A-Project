#Creates synthetic file with maximum shift of 10 million
python synthetic_data_creator.py hg19_chr2filter.bed shifted_10mtest.bed 10000000 0.1

#Find number of regions that intersect
intn=$(bedtools intersect -a hg19_chr2filter.bed -b shifted_100mtest.bed | wc -l)   

#Find number of regions that don't intersect
first=$(bedtools subtract -a hg19_chr2filter.bed -b shifted_100mtest.bed | wc -l)
second=$(bedtools subtract -a shifted_100mtest.bed -b hg19_chr2filter.bed | wc -l)

#Print ratio
echo 'Ratio:'
echo $intn/$((first + second))| bc -l
