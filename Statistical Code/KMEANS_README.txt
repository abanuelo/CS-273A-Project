Kmeans testing:
run "python scoreGenomicIntervals track1.bed track2.bed track3.bed domainfile.bed outputfile.bed bp_counting_window"
domain file is the full range over which intervals should be counted, i.e. the full length of each chromosome
run "python kmeansTest.py scoresfile.txt col1 scoresfile.txt col2 kmin kmax", optionally with "scoresfile.txt col3" at end where scoresfile.txt is output of scoreGenomicIntervals for the tracks

