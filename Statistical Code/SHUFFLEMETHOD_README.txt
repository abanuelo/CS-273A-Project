Laura's part / random shuffles and overlapSelect:

uniform_random_shuffles.py: 
	- runs all three shuffles (uniform, chromosome-restricted, gc-content-restricted) on real data
	- runs synthetic_data_creator.py on various inputs
	- runs all tests on synthetic data and plots results

make_gc_bed.py:
	- makes bed files annotated with gc-content, necessary for gc-content-restricted shuffle
	- uses fasta file, counts C,G occurences within each region