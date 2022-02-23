Each step was stored under function folder.

Step 1: Building a benchmark for motif finding
        Generate SC random sequences (with uniform frequencies of A,C,G,T) with length SL
        Generate a random motif (position weight matrix) of length ML, with total		information content being  ICPC*ML
	Generate SC string of length ML (site)
	"Plant" one sampled site at a random location in each random sequences from 		previous.
	Write out the SC sequences into a FASTA format file
	Write down the location of the planted site in each sequence in a separate text 	file called "sites,txt"
	Write down the motif that was generated previously in a separate text file called 	"motif.txt"
 	Write down the motif length in a separate test file called "motiflength.txt"

Step 2: Write a program that will read the "sequences.fa" file and "motiflength.txt" filee 	and find the motif.

Step 3: Evaluate the motif finder using KL divergence, number of overlapping positions, 	number of overlapping sites and running time.

main.py ，step1 + step2 + step3

test.py ，step2 + step3。
