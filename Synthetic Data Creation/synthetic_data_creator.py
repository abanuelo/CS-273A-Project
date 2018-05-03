
import numpy as np
import random
import sys

'''
Synthetic Data Creator: Main Contributors --> Armando Banuelos and Alex Porter

Note: original bed file comes from /afs/.ir/class/cs273a/hw1/hg18.nonexonicDeepUltras.bed 
Note 2: Right now we have symmetric alignment of bed file (ie) same value being added 
        (or substracted) from start_coord will be added (or substracted) from end_coord

To Run: python2.7 synthetic_data_creator.py hg18.nonexonicDeepUltras.bed newBedFile.bed randomshiftthreshold epislon
'''

def main():
    delimited_line = []
    threshold = 100000 #holds the max range of subtraction/addition from start/end_coord
    add_or_substract = 1 #50-50 probability --> 0 = add; 1 = subtract
    if len(sys.argv)>2:
	threshold = int(sys.argv[3])
    shiftmin = threshold/2
    epsilon = float(sys.argv[4]) #originally set to 0.1
    print "random shift range:",threshold
    print "epsilon value:", epsilon
    executeResult = 10 #make a 1/10 probability of executing nothing at all

    with open(sys.argv[1]) as fp:
        for line in fp:
            delimited_line = line.split()
            start_coord = int(delimited_line[1])
            end_coord = int(delimited_line[2])

            #Generate a random range to shift this number by
            add_or_substract = random.randint(0,1)

            #This part of the process simply returns what is still the same within
            #the original bed file
            executeResult = random.randint(1,10)
            if executeResult == 1:
                writeToNewBedFile(delimited_line)
                continue
            #Add or substract symmetrically to start_end coordinates
            if add_or_substract == 0:
                addition = random.randint(shiftmin, threshold)
                start_coord += addition
                end_coord += addition
            else:
                subtraction = random.randint(shiftmin, threshold)
                start_coord -= subtraction
                end_coord -= subtraction

            start_coord = max(start_coord,1)
            end_coord = max(end_coord,1)

            lenscale = random.uniform(-1.0*epsilon,epsilon)
            sizechange = float(lenscale)*(end_coord-start_coord)
            shiftchoice = random.randint(0,1)
        
            #only move start size if not going to hit <0
            if shiftchoice == 0 and start_coord + sizechange>0:
                start_coord = start_coord + sizechange
            else:
                end_coord = end_coord + sizechange
                start_coord = int(start_coord)
                end_coord = int(end_coord)
          
            #writes to the newBedFile.bed; 
            delimited_line[1] = int(start_coord)
            delimited_line[2] = int(end_coord)
            writeToNewBedFile(delimited_line)
            delimited_line = []

def writeToNewBedFile(delimited_line):
    #reads in file to be generated
    orig_stdout = sys.stdout
    f = open(sys.argv[2], 'a')
    sys.stdout = f

    result = "" #output to newBedFile.bed (ie chr22 2020202 3030303)
    for i in range(0, len(delimited_line)):
        if i == 0:
            result = result + str(delimited_line[i])
        else:
            result = result + "\t" + str(delimited_line[i])
    
    print result #prints result to newBedFile.bed

    #closing of the file we are outputting our results to
    sys.stdout = orig_stdout
    f.close()


if __name__ == '__main__':
    main()
