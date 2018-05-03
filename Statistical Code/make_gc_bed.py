'''
Make GC-content Annotated Bed Files
Author: Laura Miron
'''

## make bed file annotated with GC content from fasta
def GC_bed_from_fasta(fastafile='muscle_histones.fa',bedfile='muscles_gc_content.bed'):
	fasta = open(fastafile,'r').read()
	entries = fasta.split('>')[1:]
	of = open(bedfile,'w')
	for entry in entries:
		lines = entry.split('\n')
		header = lines[0]
		rang = header.split()[1]
		chrom = rang.split('=')[1].split(':')[0]
		start = rang.split(':')[1].split('-')[0]
		end = rang.split(':')[1].split('-')[1]
		seqlength = 0
		gccount = 0
		for line in lines[1:]:
			for letter in line.strip('\n'):
				seqlength += 1
				if letter in ('G','C'): gccount += 1
		percent = float(gccount)/seqlength
		of.write(chrom+'\t'+start+'\t'+end+'\t'+str(percent)+'\n')
	of.close()


# divide bed file appended with gc content into buckets, output to separate files
# (these will be input files to bedtools shuffle)
def group_by_gc_content(infile,outpattern):
	flow = open(outpattern+'_low.bed','w')
	fmid = open(outpattern+'_mid.bed','w')
	fhigh = open(outpattern+'_high.bed','w')
	with open(infile,'r') as f:
		for line in f:
			gc = float(line.split()[3]) * 100
			if gc >= 70.0:
				fhigh.write(line)
			elif gc <= 30.0:
				flow.write(line)
			else:
				fmid.write(line)
	flow.close()
	fmid.close()
	fhigh.close()


# Main Script
# make bed file of 500 bp windows of hg19, annotated with gc-content
GC_bed_from_fasta(fastafile='hg19.500_windows.fa',bedfile='gc_reference.bed')
group_by_gc_content(infile='gc_reference.bed',outpattern='gc_reference')

# make bed file of H3K27ac regions, annotated with gc-content
GC_bed_from_fasta(fastafile='muscle_histones.fa',bedfile='muscles_gc_content.bed')
group_by_gc_content(infile='muscles_gc_content.bed',outpattern='muscles_gc_content')

# make bed file of unshuffled synthetic data (hg19.chr2filter.bed)
GC_bed_from_fasta(fastafile='synthetic_gc.fa',bedfile='synthetic_gc_content.bed')
group_by_gc_content(infile='synthetic_gc_content.bed',outpattern='synthetic_gc_content')





