import numpy as np
from sklearn.cluster import KMeans
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import sys
import math
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def makeImgDir(name):
    f = name #os.path.dirname(name)
    if not os.path.exists(f):
        os.makedirs(f)

def minmaxnorm(f):
	return (f-float(f.min()))/float(f.max()-f.min())



def runKMean(X,clustercount,doprints,doplot,testpath,thirdcol):
	kmeans = KMeans(n_clusters=clustercount).fit(X)
	if doprints>0:
		print kmeans.labels_
	kmeans.cluster_centers_

	datasets = []
	averageDists = []
	centerX = []
	centerY = []
	centerZ = []
	for i in range(0,clustercount):
		datasets.append([])
		averageDists.append(0)
	for i in range(0,len(f1)):
		l = kmeans.labels_[i]
		#print "adding data set for l=",l
		if thirdcol>-1:
			datasets[l].append((f1[i],f2[i],f3[i]))
		else:
			datasets[l].append((f1[i],f2[i]))
		c1 = kmeans.cluster_centers_[l][0]
		c2 = kmeans.cluster_centers_[l][1]
		if thirdcol>-1:
			c3 = kmeans.cluster_centers_[l][2]
			centerZ.append(c3)
		centerX.append(c1)
		centerY.append(c2)
		pointdist = math.sqrt(math.fabs(c1-f1[i])**2+math.fabs(c2-f2[i])**2)
		if doprints>0:
			print "point: ",f1[i],f2[i]," cluster: ",c1,c2," dist ",pointdist	
		averageDists[l] = averageDists[l]+pointdist
		


	pointsizes =[]
	ballsizes = []
	ballmax = 0.0
	centersfile = open("centers.txt","w")
	for l in range(0,clustercount):
		if len(datasets[l])>0:
			if thirdcol>-1:
				print "center: ",centerX[l],centerY[l],centerZ[l]," size: ",len(datasets[l])
				ballsizes.append(len(datasets[l]))
				if (len(datasets[l])>ballmax):
					ballmax = len(datasets[l])
				centersfile.write(str(centerX[l])+"	"+str(centerY[l])+"	"+str(centerZ[l])+"\n")
			else:
				print "center: ",centerX[l],centerY[l]," size: ",len(datasets[l])
			averageDists[l] = averageDists[l]/float(len(datasets[l]))
			pointsizes.append(averageDists[l]*500.0)
			if doprints>0:
				print "dist sum: ",averageDists[l]," count: ",float(len(datasets[l]))
				print pointsizes[l]
		else:
			averageDists[l] = 0
			ballsizes.append(0.0)
			pointsizes.append(1.0)
	centersfile.close()
	for i in range(0,len(ballsizes)):
		ballsizes[i] = float(ballsizes[i])/ballmax*2.0	
		if ballsizes[i] ==0:
			ballsizes[i] = 0.1
	print ballsizes
	score = np.mean(np.array(averageDists))
	lseries = pd.Series((kmeans.labels_[i] for i in range(0,len(f1))))
	df2 = df.assign(l=lseries)
	do3d = True
	if doplot>0:
		if thirdcol>-1:
			if do3d:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				ax.scatter(centerZ,centerY,centerX,c = 'r',s= ballsizes)
				#ax.scatter(f3,f2,f1,c = df2['l'],s=0.5)
				ax.set_zlabel('Species Conservation')
				ax.set_ylabel('H3K27ac Conserved Modifications')
				ax.set_xlabel('H3K27me3 Conserved Modifications')

				#plt.show()
				for angle in range(0, 360,30):
				    ax.view_init(30, angle)
				    plt.draw()
				    plt.savefig(testpath+"_"+str(clustercount)+"_"+str(angle)+"center.png")
				for angle in range(0, 360,30):
				    ax.view_init(angle, 30)
				    plt.draw()
				    plt.savefig(testpath+"_"+str(clustercount)+"_"+str(angle)+"bcenter.png")
						   
				plt.close()
			else:
				plt.scatter(f2,f3,c=f1,s= 0.01)
				plt.savefig(testpath+"_"+str(clustercount)+".png")
				plt.close()
		else:
			plt.scatter(f1,f2,c= df2['l'])
			plt.scatter(centerX,centerY,c='r')
			plt.xlabel('Species Conservation')
			#plt.ylabel('H3K27ac Conserved Modifications')
			plt.ylabel('H3K27me3 Conserved Modifications')
			plt.savefig(testpath+"_"+str(clustercount)+".png")
			plt.close()
	return score


args = sys.argv
if len(sys.argv)!=7 and len(sys.argv)!=8:
	print "ERR: format is script.py datafile col1 datafile2 col2 kmin kmax or {col2} ->{difcol1 difcol2}"


chromcolnum = 0
filename1 = args[1]
filename2 = args[3]
kmin = int(args[5])
kmax = int(args[6])
col1 = int(args[2])
col2 = int(args[4]) 
col3 = -1
filename3 = ""
if len(sys.argv)>8:
	filename3 = args[7]
	col3 = int(args[8])
col2m = -1

df = pd.read_csv(filename1, sep='\t',header=None)
print df
df2 = pd.read_csv(filename2,sep = '\t',header = None)
#df2[(np.abs(stats.zscore(df2)) < 3).all(axis=1)]
#df[(np.abs(stats.zscore(df)) < 3).all(axis=1)]

print df2
df3 = df2
if col3>-1:
	df3 = pd.read_csv(filename3,sep='\t',header = None)

f1 =  df[col1].values
f2 =  df2[col2].values
f3 = f2
if col3>-1:
	f3 = df3[col3].values
f1 = f1+0.01
f1[0] = f1[0]+0.02
f2 = f2+1
print f1
print f2
print f3
chromcols = df[chromcolnum].values
#print chromcols
chromnumeric = chromcols
if chromcols[0].find("chr")>0:
	chromnumeric = []
	for name in chromcols:
		tempname = name.replace("chr","")
		tempname = tempname.replace("X","23")
		tempname = tempname.replace("Y","24")
		chromnumeric.append(int(tempname))
	
if col1 == chromcolnum:
	f1 = np.array(chromnumeric)
if col2 == chromcolnum:
	f2 = np.array(chromnumeric) 

f1 = minmaxnorm(f1)
f2 = minmaxnorm(f2)
if col3>-1:
	f3 = minmaxnorm(f3)
X = np.matrix(zip(f1,f2))
if col3>-1:
	X = np.matrix(zip(f1,f2,f3))
scores = []
foldername = "kmeantests"
makeImgDir(foldername)
for k in range(kmin,kmax):
	score = runKMean(X,k,0,1,foldername+"/clusterplot",col3)
	print "score: ",score
	scores.append(score)


#plt.figure()
plt.plot(range(kmin,kmax),scores)
plt.scatter(range(kmin,kmax),scores)
plt.show()
