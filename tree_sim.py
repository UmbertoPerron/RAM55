#!/usr/bin/env python	
from __future__ import with_statement

from sys import stdin, stderr, stdout
from optparse import OptionParser
from vfork.io.util import safe_rstrip
from vfork.util import exit, format_usage
#from random import shuffle
#from subprocess import Popen, PIPE
#from collections import defaultdict
#from vfork.io.colreader import Reader

from ete2 import Tree
import numpy as np
import distance
import scipy.linalg
import math
import random

import pyvolve

masking_dict = {'1': 'V', '0': 'V', '2': 'V', 'A': 'A', 'C': 'R', 'B': 'R', 'E': 'N', 'D': 'R', 'G': 'N', 'F': 'N', 'I': 'D', 'H': 'D', 'K': 'C', 'J': 'D', 'M': 'C', 'L': 'C', 'O': 'Q', 'N': 'Q', 'Q': 'E', 'P': 'Q', 'S': 'E', 'R': 'E', 'U': 'H', 'T': 'G', 'W': 'H', 'V': 'H', 'Y': 'I', 'X': 'I', 'Z': 'I', 'a': 'L', 'c': 'L', 'b': 'L', 'e': 'K', 'd': 'K', 'g': 'M', 'f': 'K', 'i': 'M', 'h': 'M', 'k': 'F', 'j': 'F', 'm': 'P', 'l': 'F', 'o': 'S', 'n': 'P', 'q': 'S', 'p': 'S', 's': 'T', 'r': 'T', 'u': 'W', 't': 'T', 'w': 'W', 'v': 'W', 'y': 'Y', 'x': 'Y', 'z': 'Y'}
		

def main():
	usage = format_usage('''
		%prog NEWICK_FILE LABELS < ORIG_MATRIX 
	''')
	parser = OptionParser(usage=usage)
	
	parser.add_option('-m', '--mask', default=False, action='store_true', help='mask a rotasequence to make it look like a AA sequnce (1-character symbols) [default: %default]')
	parser.add_option('-u', '--uniform_model', default=False, action='store_true', help='simulate sequence evolution under a uniform model [default: %default]')
	parser.add_option('-f', '--file_seq', default=False, help='custom sequence, symbols must be space-separated [default: %default]')
	parser.add_option('-r', '--runs', type="int", default=False, help='number of simulation runs to be performed [default: 5]')
	parser.add_option('-F', '--Freq_file', default=False, help='use frequencies from file to generate the original (i.e. ancestor) sequence [default: %default]')
	parser.add_option('-l', '--length', type="int", default=False, help='length of sequence to be generated [default length: 1000]')
	parser.add_option('-g', '--gamma', type="float", default=False, help='use heterogeneous rates drawn from a gamma distribution of shape a (with a = b and mean=1)')
	parser.add_option('-o', '--outfile',  default=stdout, help='output filename pattern for multiple runs in separate files [default length: 1000]')
	parser.add_option('-p', '--pyvolve',  default=False, action='store_true', help='use pyvolve to simulate evolution [default: %default]')
	options, args = parser.parse_args()
	
	if len(args) != 2:
		exit('Unexpected argument number.')
	elif options.pyvolve and not options.Freq_file:
		exit('pyvolve requires the equilibrium frequencies')
	# can directly generate sequences using 1-character alphabet
	LABELS = args[1].split(",")
	nstates = len(LABELS)

	orig_matrix = np.zeros((nstates,nstates), dtype='float64')
	if options.uniform_model:
		# generate a uniform model to simulate from
		for i in range(nstates):
			# all rows in rotamers_adj_count_matrix.merged.Dayhoff_norm.IRM.Q-star sum to 5.10078 
			# 5.10078 / 55 = 0.09274
			orig_matrix[i,:] = [0.09274] * 55
	else:
		# parse the original IRM matrix
		i = 0
		for line in stdin:
			orig_matrix[i,:] = safe_rstrip(line).split('\t')
			i = i + 1
	if not options.pyvolve:
		# from IRM to Q with lines that sum to 0
		for i in range(0,nstates):
			orig_matrix[i,i] = (np.sum(orig_matrix[i,:]) - orig_matrix[i,i])*-1
		tree = Tree(args[0], format=1)

	# specify sequence length if one has to be generated
	if options.length:
		seq_len = options.length
	else:
		seq_len = 1000

	# either read the ancestor sequence from file
	if options.file_seq:
		with open(options.file_seq, 'r') as fs:
			for line in fs:
				orig_seq = safe_rstrip(line).split(' ')
		seq_len = len(orig_seq)
		for symbol in orig_seq:
			if symbol not in LABELS:
				exit('Unexpected symbol %s found in file_seq') % (symbol)
	# or generate it according to specified frequencies
	if options.Freq_file:
		with open(options.Freq_file, 'r') as fs:
			# remove header
			lines = fs.readlines()
			freqs = [float(x) for x in safe_rstrip(lines[1]).split('\t')]
			if len(freqs) != len(LABELS):
				exit('Submitted frequencies should match states in LABELS')
		if not options.file_seq:
			orig_seq = np.random.choice(LABELS, seq_len, p=freqs)
	# or generate a truly random sequence
	elif not options.Freq_file and not options.file_seq:
		orig_seq = np.random.choice(LABELS, seq_len)


	# specify the number of simulation runs i.e the number of alignments to be generated 
	# starting from the same ancestor sequence
	if options.runs:
		runs = options.runs
	else:
		runs = 1

	if options.pyvolve:
		# Define a custom Q matrix "This rate matrix must be square and all rows in this matrix must sum to 0."
		# https://github.com/sjspielman/pyvolve/blob/master/user_manual/pyvolve_manual.pdf
		custom_matrix = orig_matrix 

		# Define custom frequencies
		#fi = pyvolve.CustomFrequencies("", freq_dict = dict(zip(LABELS, freqs))) 
		#frequencies = f.compute_frequencies()
		#frequencies = freqs

		# Construct a model using the custom rate-matrix,code and frequencies
		#custom_model = pyvolve.Model("custom", {"matrix":custom_matrix, "code":LABELS, "state_freqs":frequencies})	
		custom_model = pyvolve.Model("custom", {"matrix":custom_matrix, "code":LABELS})

		# Define a Partition object with a specified ancestral sequence
		my_partition = pyvolve.Partition(models = custom_model, root_sequence = ''.join(map(str, orig_seq)))

		# Define an Evolver instance to evolve a single partition
		tree = pyvolve.read_tree(file = args[0])
		my_evolver = pyvolve.Evolver(partitions = my_partition, tree = tree) 
	
	if options.gamma:
		# uniquely assign a rate r to every site drawing from a gamma distribution
		# with shape = alpha, scale( aka theta ) = 1 / beta = 1 / alpha, and mean = 1
		shape = float(options.gamma)
		scale = 1 / shape
		gamma_rates = np.random.gamma(shape, scale, seq_len)

	for run_num in range(runs):
		if options.outfile != stdout:
			outfile = open(options.outfile + '_' + str(run_num), "w")
		else:
			outfile = stdout
		if options.pyvolve:
			my_evolver(seqfile = outfile, seqfmt = "phylip-sequential", ratefile = None, infofile = None)
		else:
			# generate a dictionary of sequences accessed by node names
			seq_dict = {}
			# assign the ancestor sequence to the root node 
			seq_dict[tree.get_tree_root().name] = orig_seq

			# traverse the tree using a preorder strategy (avoiding the root node):
			# 1) Visit the first descendant node, 2) Traverse the left subtree , 3) Traverse the right subtree
			for node in tree.iter_descendants("preorder"):
				# calculate P(t), evolve the parent sequence accordingly
				parent_seq = seq_dict[node.up.name]
				t = node.dist
				new_seq = np.zeros(seq_len, dtype=object)
				if not options.gamma:
					Pt_matrix = scipy.linalg.expm(t* orig_matrix)
				for site,current_state in enumerate(parent_seq):
					if options.gamma:
						t_r = t * gamma_rates[site]
						Pt_matrix = scipy.linalg.expm(t_r* orig_matrix)
					P_state = Pt_matrix[LABELS.index(current_state),:]
					new_state = np.random.choice(LABELS, 1, p=P_state)[0]
					new_seq[site] = new_state
				seq_dict[node.name] = new_seq
			
			# iterate over the leaves, print a PHYLIP alignment
			leaves = [leaf.name for leaf in tree.iter_leaves()]
			print >>outfile, "%d\t%d" % (len(leaves), len(orig_seq))
			print >>outfile, ''
			for leaf in leaves:
				sequence = seq_dict[leaf]

				# mask the rotamer states
				if options.mask:
					masked_sequence = np.zeros(seq_len, dtype=object)
					for site, state in sequence:
						masked_sequence[site] = masking_dict[state]
					sequence = masked_sequence
				
				print >>outfile, "%s\t%s" % (leaf, ' '.join(map(str, sequence))) 	

if __name__ == '__main__':
	main()

