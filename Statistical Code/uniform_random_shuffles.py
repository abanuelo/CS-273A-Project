'''
Uniform Random Shuffles
Author: Laura Miron
'''

import os
import numpy as np
import pdb
import matplotlib.pyplot as plt
import pickle

histone_file = 'muscles.short.bed'
conserved_file = 'consindelsHgMmCanFam.bed'
num_shuffles = 10
true_overlapsfile = 'trueOverlaps.bed'
true_countsfile = 'counts.txt'

# run overlapSelect on specified files and return number of overlaps
def overlapSelect(select_file,infile,outfile):
	os.system(r'/afs/ir/class/cs273a/bin/overlapSelect_64 '+select_file+' '+infile+' '+outfile)
	return len(open(outfile,'r').readlines())

# perform uniform bedtools shuffle on specified files for num_shuffles iterations
# on each iteration, perform overlapSelect on input_file and shuffled file and output number of overlaps to shuffled_countsfile
def uniform_shuffle(select_file=conserved_file,input_file=histone_file,shuffled_countsfile='shuffledCounts.txt'):
	print 'starting uniform shuffle'
	if os.path.isfile(shuffled_countsfile):
		os.remove(shuffled_countsfile)
	for i in range(0,num_shuffles):
		shufflefile = input_file+'.uniform.'+str(i)+'.bed'
		if os.path.isfile(shufflefile):
			os.remove(shufflefile)
		os.system(r'bedtools shuffle -i '+input_file+' -g chrom.sizes | cut -f1-3 > '+shufflefile)

		overlapfile = 'uniformOverlaps.'+str(i)+'.bed'
		num_overlaps = overlapSelect(select_file=select_file,infile=shufflefile,outfile=overlapfile)
		os.system(r'echo ' +str(num_overlaps)+' >> '+shuffled_countsfile)

# perform bedtools shuffle restricted by chromosome on specified files for num_shuffles iterations
# on each iteration, perform overlapSelect on input_file and shuffled file and output number of overlaps to shuffled_countsfile
def chromosome_shuffle(select_file=conserved_file,input_file=histone_file,shuffled_countsfile='chromosomeShuffledCounts.txt'):
	print 'starting chromosome-restricted shuffle'
	if os.path.isfile(shuffled_countsfile):
		os.remove(shuffled_countsfile)
	for i in range(0,num_shuffles):
		shufflefile = input_file+'.chromosome.'+str(i)+'.bed'
		if os.path.isfile(shufflefile):
			os.remove(shufflefile)
		os.system(r'bedtools shuffle -i '+input_file+' -g chrom.sizes | cut -f1-3 > '+shufflefile)

		overlapfile = 'chromosomeOverlaps.'+str(i)+'.bed'
		num_overlaps = overlapSelect(select_file=select_file,infile=shufflefile,outfile=overlapfile)
		os.system(r'echo ' +str(num_overlaps)+' >> '+shuffled_countsfile)

# perform bedtools shuffle restricted by gc_content on specified files for num_shuffles iterations
# on each iteration, perform overlapSelect on input_file and shuffled file and output number of overlaps to shuffled_countsfile
def gc_content_shuffle(select_file=conserved_file,inpattern='muscles_gc_content',shuffled_countsfile='gcShuffledCounts.txt'):
	print 'starting gc content-restricted shuffle'
	if os.path.isfile(shuffled_countsfile): os.remove(shuffled_countsfile)
	for i in range(0,num_shuffles):
		lowshufflefile = 'shuffled_gc_low.'+str(i)+'.bed'
		midshufflefile = 'shuffled_gc_mid.'+str(i)+'.bed'
		highshufflefile = 'shuffled_gc_high.'+str(i)+'.bed'
		fullshufflefile = 'shuffled_gc_full.'+str(i)+'.bed'

		# -incl files are bed files of 500 bp windows across entire hg19 genome, sorted by gc percentage into three files
		os.system(r'bedtools shuffle -i '+inpattern+'_low.bed -g chrom.sizes -incl gc_reference_low.bed > '+lowshufflefile)
		os.system(r'bedtools shuffle -i '+inpattern+'_mid.bed -g chrom.sizes -incl gc_reference_mid.bed > '+midshufflefile)
		os.system(r'bedtools shuffle -i '+inpattern+'_high.bed -g chrom.sizes -incl gc_reference_high.bed > '+highshufflefile)
		os.system(r'cat '+lowshufflefile+' '+midshufflefile+' '+highshufflefile+' > '+fullshufflefile)
		
		overlapfile = 'gcOverlaps.'+str(i)+'.bed'
		num_overlaps = overlapSelect(select_file=select_file,infile=fullshufflefile,outfile=overlapfile)
		os.system(r'echo ' +str(num_overlaps)+' >> '+shuffled_countsfile)

# calculate fold and zscore
def calculateFoldAndZscore(true_counts,shuffled_counts):
	shuffled_avg = float(shuffled_counts.sum()) / len(shuffled_counts)
	std_dev = np.std(shuffled_counts)
	fold = float(true_counts) / shuffled_avg
	zscore = float(true_counts - shuffled_avg) / std_dev
	return (fold,zscore)

# read all overlap counts from shuffled_countsfile and true_countsfile, calculate fold and zscore
def printFoldAndZscore(shuffletype,shuffled_countsfile,true_countsfile=true_countsfile):
	print 'calculating fold and zscore'
	true_counts = int(open(true_countsfile).readlines()[0].rstrip('\n'))
	shuffled_counts = np.array([int(line.rstrip('\n')) for line in open(shuffled_countsfile)])
	fold,zscore = calculateFoldAndZscore(true_counts=true_counts,shuffled_counts=shuffled_counts)
	print shuffletype + " fold:" + str(fold) + ", zscore:" + str(zscore) + "\n"

# invoke synthetic data creator with different thresholds
def make_synthetic_data():
	for threshold in range(100,1000,100):
		epsilon = 0.1
		outfile = 'hg19chr2_'+str(threshold)+'_p'+str(epsilon)[2:]+'.bed'
		cmd = 'python synthetic_data_creator.py hg19_chr2filter.bed '+outfile+' '+str(threshold)+' '+str(epsilon)
		print cmd
		os.system(cmd)

# run all tests (uniform_shuffle, chromosome_shuffle, gc_content_shuffle) on synthetic data	files
# in range specified by threshold, output results to files to be plotted	
def test_synthetic_data():
	print 'testing synthetic data'
	infile = 'hg19_chr2filter.bed'
	xaxis = []
	uniform_fold = []
	uniform_zscore = []
	chromosome_fold = []
	chromosome_zscore = []
	gc_content_fold = []
	gc_content_zscore = []
	
	for threshold in range(10000,100000,10000): # only makes sense if make_synethetic data have previous been run with these values
		xaxis.append(threshold)
		epsilon = 0.1
		synfile = 'hg19chr2_'+str(threshold)+'_p'+str(epsilon)[2:]+'.bed'
		true_counts = overlapSelect(select_file=infile,infile=synfile,outfile='syntheticTrueCounts.txt')
		
		# perform uniform shuffle on synthetic data
		uniform_shuffled_countsfile='shuffledSyntheticUniform.txt'
		uniform_shuffle(select_file=infile,input_file=synfile,shuffled_countsfile=uniform_shuffled_countsfile)
		shuffled_counts = np.array([int(line.rstrip('\n')) for line in open(uniform_shuffled_countsfile)])
		fold,zscore = calculateFoldAndZscore(true_counts=true_counts,shuffled_counts=shuffled_counts)
		uniform_fold.append(fold)
		uniform_zscore.append(zscore)
		print 'uniform synthetic '+str(threshold)+' '+str(epsilon)+' fold: '+str(fold)+' zscore: '+str(zscore)

		# perform chromosome-restricted shuffle on synthetic data
		chromosome_shuffled_countsfile='shuffledSyntheticChromosome.txt'
		chromosome_shuffle(input_file=infile,select_file=synfile,shuffled_countsfile=chromosome_shuffled_countsfile)
		shuffled_counts = np.array([int(line.rstrip('\n')) for line in open(chromosome_shuffled_countsfile)])
		fold,zscore = calculateFoldAndZscore(true_counts=true_counts,shuffled_counts=shuffled_counts)
		chromosome_fold.append(fold)
		chromosome_zscore.append(zscore)
		print 'chromosome synthetic '+str(threshold)+' '+str(epsilon)+' fold: '+str(fold)+' zscore: '+str(zscore)

		# perform gc-content-restricted shuffle on synthetic data
		gc_content_shuffled_countsfile='shuffledSyntheticChromosome.txt'
		gc_content_shuffle(select_file=synfile,inpattern='synthetic_gc_content',shuffled_countsfile=gc_content_shuffled_countsfile)
		shuffled_counts = np.array([int(line.rstrip('\n')) for line in open(gc_content_shuffled_countsfile)])
		fold,zscore = calculateFoldAndZscore(true_counts=true_counts,shuffled_counts=shuffled_counts)
		gc_content_fold.append(fold)
		gc_content_zscore.append(zscore)
		print 'gc_content synthetic '+str(threshold)+' '+str(epsilon)+' fold: '+str(fold)+' zscore: '+str(zscore)

	# output results to files
	with open('xaxis.txt','w') as f: pickle.dump(xaxis,f)
	with open('uniform_fold.txt','w') as f: pickle.dump(uniform_fold,f)
	with open('uniform_zscore.txt','w') as f: pickle.dump(uniform_zscore,f)
	with open('chromosome_fold.txt','w') as f: pickle.dump(chromosome_fold,f)
	with open('chromosome_zscore.txt','w') as f: pickle.dump(chromosome_zscore,f)
	with open('gc_content_fold.txt','w') as f: pickle.dump(gc_content_fold,f)
	with open('gc_content_zscore.txt','w') as f: pickle.dump(gc_content_zscore,f)

# make plots of fold and zscore
def plot_synthetic_data():
	# read data output by test_synthetic_data
	with open('xaxis.txt','rb') as f: xaxis = pickle.load(f)
	with open('uniform_fold.txt','rb') as f: uniform_fold = pickle.load(f)
	with open('chromosome_fold.txt','rb') as f: chromosome_fold = pickle.load(f)
	with open('gc_content_fold.txt','rb') as f: gc_content_fold = pickle.load(f)
	with open('uniform_zscore.txt','rb') as f: uniform_zscore = pickle.load(f)
	with open('chromosome_zscore.txt','rb') as f: chromosome_zscore = pickle.load(f)
	with open('gc_content_zscore.txt','rb') as f: gc_content_zscore = pickle.load(f)

	# plot fold vs avg threshold
	xaxis = [int(x)/2 for x in xaxis]
	plt.plot(xaxis,uniform_fold,label='uniform shuffle')
	plt.plot(xaxis,chromosome_fold,label='chromosome-restricted shuffle')
	plt.plot(xaxis,gc_content_fold,label='gc-content-restricted shuffle')
	plt.legend()
	plt.title('Fold Reported by Shuffles on Synthetic Data')
	plt.xlabel('avg shift per bed line (bp)')
	plt.ylabel('fold')
	plt.show()

	# plot zscore vs avg threshold
	plt.plot(xaxis,uniform_zscore,label='uniform shuffle')
	plt.plot(xaxis,chromosome_zscore,label='chromosome-restricted shuffle')
	plt.plot(xaxis,gc_content_zscore,label='gc-content-restricted shuffle')
	plt.title('Zscore Reported by Shuffles on Synthetic Data')
	plt.xlabel('avg shift per bed line (bp)')
	plt.ylabel('zscore')
	plt.show()


# Main Script

# make and test synthetic data
make_synthetic_data()
test_synthetic_data()
plot_synthetic_data()

# run overlap select on unshuffled H3K27ac file and species-conserved regions file
os.system(r'/afs/ir/class/cs273a/bin/overlapSelect_64 '+conserved_file+' '+histone_file+' '+true_overlapsfile)
os.system(r'cat ' + true_overlapsfile + ' | wc -l > '+true_countsfile)

# perform uniform shuffle, chromosome-restricted shuffle, and gc-content-restricted shuffle on 
# H3K27ac and species-conserved regions, output fold and zscore
uniform_shuffle()
printFoldAndZscore(shuffletype='uniform',shuffled_countsfile='shuffledCounts.txt')
chromosome_shuffle()
printFoldAndZscore(shuffletype='chromosome-restricted',shuffled_countsfile='chromosomeShuffledCounts.txt')
gc_content_shuffle()
printFoldAndZscore(shuffletype='gc-content-restricted',shuffled_countsfile='gcShuffledCounts.txt')



