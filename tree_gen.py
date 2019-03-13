#!/usr/bin/env python	
from __future__ import with_statement

from sys import stdin, stderr
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
import math
import random
import string

def main():
	usage = format_usage('''
		%prog LEAVES_NUM
		
		Generate a random newick tree
	''')
	parser = OptionParser(usage=usage)

	parser.add_option('-b', '--branch_range', default=False, help='custom branch length range, values should be comma separated [default range: 0.01,2]')
	options, args = parser.parse_args()
	
	if len(args) != 1:
		exit('Unexpected argument number.')

	leaves_num = int(args[0])
	nodes_num = leaves_num*2 
	symbols = string.ascii_letters
	# a list of node names
	nodes_names = [''.join([random.choice(symbols) for i in range(5)]) for j in range(nodes_num)]
	if options.branch_range:
		branch_range = tuple([float(x) for x in options.branch_range.split(",")])
	else:
		branch_range = (0.01, 2)

	tree = Tree()
	# generate the tree
	tree.populate(leaves_num, random_branches=True, branch_range=branch_range)
	
	# name all nodes, the populate names_library function only names branches
	i = 0
	for node in tree.traverse("preorder"):
		node.name = nodes_names[i]
		i = i+1

	# print it in newick format, include the internal node names
	print tree.write(format=1)

if __name__ == '__main__':
	main()

