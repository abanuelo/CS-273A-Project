#create bed file with scores for conservation, and some histone mods to kmeans-cluster
#for line in conservation file, put conservation score and # of each histone type overlaps
#use: script.py conservationfile modfile1 modfile2 chromosomerangesfile testinginterval



def mapBins(filename,baseinterval,chrset):
	consmap = {}
	for chrname in chrset:
		consmap[chrname] = {}
	with open(filename,"r") as f:
		for line in f:
			dat = line.split()
			start = int(dat[1])
			chrname = dat[0]
			startbin = start - (start % baseinterval)
			startid = startbin
			#print start,startbin
			if startid not in consmap[chrname].keys():
				consmap[chrname][startid] = 0
			consmap[chrname][startid]+=1

	return consmap






import sys
import gzip
from difflib import SequenceMatcher
args = sys.argv
if len(sys.argv)!=3:
        print "ERR: format is script.py infile outfile"
consfile = str(args[1])
mod1 = str(args[2])
mod2 = str(args[3])
chrrangefile = str(args[4])
outfile = str(args[5])
baseinterval = int(args[6])
#first count # lines for each chr start end in each intersect file

chrmap = {}
chrset = []
with open(chrrangefile,"r") as f:
	for line in f:
		dat = line.split()
		chrname = dat[0]
		chrsize = int(dat[2])	
		chrset.append(chrname)
		chrmap[chrname] =chrsize		

map1 = mapBins(consfile,baseinterval,chrset)
map2 = mapBins(mod1,baseinterval,chrset)
map3 = mapBins(mod2,baseinterval,chrset)
overlaps = 0
ones = 0
twos = 0
empties = 0
outbed = open(outfile,"w+")
for chrname in chrset:
	for i in range(0,int(float(chrmap[chrname])/float(baseinterval))):
		bin = i*baseinterval
		score1 = 0
		score2= 0
		score3 = 0
		if bin in map1[chrname].keys():
			score1 = map1[chrname][bin]
		if bin in map2[chrname].keys():
			score2 = map2[chrname][bin]		
		if bin in map3[chrname].keys():
			score3 = map3[chrname][bin]
		#print score1,score2,score3
		#if (score1!=0 or score2!=0 or score3!=0) and (score1< 1000) and (score2<1000):
		if score1==0 and score2==0:
			empties =empties +1
		if score1==0 and score2!=0:
			ones = ones  +1
		if score2==0 and score1!=0:
			twos = twos+1
		if score1!=0 and score2!=0:
			overlaps = overlaps+2		
		if (score2  <100) and (score1<1000):
			outbed.write(chrname)
			outbed.write("	")
			outbed.write(str(bin))
			outbed.write("	")
			outbed.write(str(score1))
			outbed.write("	")
			outbed.write(str(score2))
			outbed.write("	")
			outbed.write(str(score3))
			outbed.write("\n")		

print "overlap: ",overlaps
print "ones: ",ones
print "twos: ",twos
print "empties: ",empties
