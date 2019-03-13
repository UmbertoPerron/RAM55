#!/usr/bin/env python
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
from vfork.io.util import safe_rstrip
from vfork.util import exit, format_usage

from collections import namedtuple
import numpy as np
import scipy.linalg
import pandas as pd
from dendropy import Tree, TaxonNamespace, Taxon
import random
import string
import warnings
import itertools

# a custom NODE obj
NODE = namedtuple('NODE', ['Cx_dict', 'Lx_dict'])
#Lx(i)is the likelihood of the best reconstruction of this subtree 
#on the condition that the father
#of node x (the current node) is assigned character state i.
#Cx(i) is the character state assigned to node x (the current node) 
#in this optimal conditional reconstruction.

one_rota_to_four_dict = {'1': 'VAL2', '0': 'VAL1', '2': 'VAL3', 'A': 'ALA', 'C': 'ARG2', 'B': 'ARG1', 'E': 'ASN1', 'D': 'ARG3', 'G': 'ASN3', 'F': 'ASN2', 'I': 'ASP2', 'H': 'ASP1', 'K': 'CYS1', 'J': 'ASP3', 'M': 'CYS3', 'L': 'CYS2', 'O': 'GLN2', 'N': 'GLN1', 'Q': 'GLU1', 'P': 'GLN3', 'S': 'GLU3', 'R': 'GLU2', 'U': 'HIS1', 'T': 'GLY', 'W': 'HIS3', 'V': 'HIS2', 'Y': 'ILE2', 'X': 'ILE1', 'Z': 'ILE3', 'a': 'LEU1', 'c': 'LEU3', 'b': 'LEU2', 'e': 'LYS2', 'd': 'LYS1', 'g': 'MET1', 'f': 'LYS3', 'i': 'MET3', 'h': 'MET2', 'k': 'PHE2', 'j': 'PHE1', 'm': 'PRO1', 'l': 'PHE3', 'o': 'SER1', 'n': 'PRO2', 'q': 'SER3', 'p': 'SER2', 's': 'THR2', 'r': 'THR1', 'u': 'TRP1', 't': 'THR3', 'w': 'TRP3', 'v': 'TRP2', 'y': 'TYR2', 'x': 'TYR1', 'z': 'TYR3'}

four_rota_to_one_dict = dict([[v,k] for k,v in one_rota_to_four_dict.items()])

three_AA_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'TRP': 'W', 'PRO': 'P', 'THR': 'T', 'ILE': 'I', 'ALA': 'A', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'MET': 'M', 'GLU': 'E', 'ASN': 'N', 'TYR': 'Y', 'VAL': 'V'}

one_AA_to_three_dict = dict([[v,k] for k,v in three_AA_to_one.items()])

three_AA_to_three_rota_dict = {'CYS': ['CYS2', 'CYS1', 'CYS3'], 'ASP': ['ASP1', 'ASP3', 'ASP2'], 'SER': ['SER2', 'SER1', 'SER3'], 'GLN': ['GLN1', 'GLN3', 'GLN2'], 'LYS': ['LYS1', 'LYS3', 'LYS2'], 'ILE': ['ILE1', 'ILE3', 'ILE2'], 'PRO': ['PRO2', 'PRO1'], 'THR': ['THR1', 'THR3', 'THR2'], 'PHE': ['PHE1', 'PHE3', 'PHE2'], 'ALA': ['ALA'], 'GLY': ['GLY'], 'HIS': ['HIS2', 'HIS1', 'HIS3'], 'GLU': ['GLU2', 'GLU1', 'GLU3'], 'LEU': ['LEU2', 'LEU1', 'LEU3'], 'ARG': ['ARG1', 'ARG3', 'ARG2'], 'TRP': ['TRP2', 'TRP1', 'TRP3'], 'VAL': ['VAL1', 'VAL3', 'VAL2'], 'ASN': ['ASN2', 'ASN1', 'ASN3'], 'TYR': ['TYR1', 'TYR3', 'TYR2'], 'MET': ['MET2', 'MET1', 'MET3']}

def AA_2_rota(one_char_AA_state):
	rota_states = three_AA_to_three_rota_dict[one_AA_to_three_dict[one_char_AA_state]]
	return [four_rota_to_one_dict[r] for r in rota_states]

def main():
	usage = format_usage('''
		%prog NEWICK_TREE PHYLIP_ALIGN IRM_SCALED FREQ_FILE

		reconstruct internal node sequences using
		the Joint Reconstruction algorithm from Pupko et al. 2000
		http://www.ncbi.nlm.nih.gov/pubmed/10833195
		''')
	parser = OptionParser(usage=usage)
	
	parser.add_option('-s', '--state_set', type=str, help='specify state set if IRM or freqs labels  and alignment symbols dont match')
	parser.add_option('-e', '--exhaustive', action="store_true", default=False, help='Iter over all states not only those observed at site')
	parser.add_option('-c', '--constrained', default=False, help='Constrain rotamer states reconstruction using internal AA sequences.')
	parser.add_option('-i', '--ignore_leaf', type=str, default=False, help='Ignore leaf sequence (i.e. treat leaf as internal node)')
	parser.add_option('-I', '--Ignore_and_trim_sibling', type=str, default=False, help='Ignore leaf sequence and do not use Lx from its sibling (i.e. treat leaf as root)')
	options, args = parser.parse_args()
	
	if len(args) != 4:
		exit('Unexpected argument number.')
	
	treefile, alignfile, IRMfile, freqfile = args

	if options.Ignore_and_trim_sibling:
		options.ignore_leaf = options.Ignore_and_trim_sibling
	
	tns = TaxonNamespace(is_case_sensitive=True)

	# dedropy creates a fake root if tree is unrooted
	tree = Tree.get(path=treefile, schema='newick', taxon_namespace=tns, case_sensitive_taxon_labels=True, suppress_internal_node_taxa=False)
	
	# assign taxon instances (and corresponding randomly-generated labels) 
	# to internal nodes if not present
	taxonless_intnodes = [intnode for intnode in tree.postorder_internal_node_iter() if intnode.taxon is None and intnode._parent_node is not None]
	if len(taxonless_intnodes) >0:
		warnings.warn('Internal node taxon labels not found in %s, generating random labels...' % (treefile))		
		symbols = string.ascii_letters
		intnode_names = [''.join([random.choice(symbols) for i in range(5)]) for j in range(len(taxonless_intnodes))]
		for intnode, intnode_name in zip(taxonless_intnodes, intnode_names):
			intnode.taxon = Taxon(label=intnode_name)
			warnings.warn('Assigned label %s to node %s' % (intnode.taxon.label, str(intnode)))
		root = [intnode for intnode in tree.postorder_internal_node_iter() if intnode._parent_node is None][0]
		root.taxon = Taxon(label='ROOT')
		warnings.warn('Assigned label %s to node %s' % (root.taxon.label, str(root)))
	internal_taxons = ['ROOT'] + [intnode.taxon.label for intnode in tree.postorder_internal_node_iter() if intnode._parent_node is not None]	

	# read alignment
	orig_align_df = pd.read_csv(alignfile,  sep="\t", header=None, index_col=0, names=["sequence"], skiprows=1)
	taxa_num, align_len = [len(orig_align_df), orig_align_df["sequence"].str.split(' ').str.len().values[0]]
	# split sequences into alignment columns, one per site
	split_align_df = orig_align_df['sequence'].str.split(' ', expand=True).rename(columns = lambda x: str(x+1))

	if options.constrained:
		# read the internal AA sequence alignment
		constraint_align_df = pd.read_csv(options.constrained,  sep="\t", header=None, index_col=0, names=["sequence"], skiprows=1)
		constraint_split_align_df = constraint_align_df['sequence'].str.split(' ', expand=True).rename(columns = lambda x: str(x+1))
	
	if options.ignore_leaf:
		reconstructed_split_align_df = pd.DataFrame(0, index=internal_taxons + [options.ignore_leaf], columns=[])
	else:
		# output alignment of reconstructed sequences from internal nodes
		reconstructed_split_align_df = pd.DataFrame(0, index=internal_taxons, columns=[]) 
	
				
	# read IRM, needs rows that sums to 0, also needs to be scaled 
	# so that it has on avg 1 AA state change per unit of time at equilibrium

	# TODO add option to process exchangeabilities instead
	if options.state_set:
		# custom state labels
		all_states = options.state_set.split(",")
		IRM_df = pd.read_csv(IRMfile, sep="\t", header=None, skiprows=1, names=all_states, index_col=0)
		IRM_df.index = all_states
		freqs_df = pd.read_csv(freqfile, sep="\t", header=None, skiprows=1, names=all_states, index_col=None)
	else:
		IRM_df = pd.read_csv(IRMfile, sep="\t", header=0, index_col=0)
		all_states = IRM_df.index.tolist()
		freqs_df = pd.read_csv(freqfile, sep="\t", header=0, index_col=None)
	nstates = len(all_states)

	# from scaled IRM to Q with lines that sum to 0
	#exponentiate Q for every brlen in tree
	# head_node is the child node for this edge
	Pt_dict = {}
	IRM_arr = IRM_df.values
	for i in range(0,nstates):
		IRM_arr[i,i] = (np.sum(IRM_arr[i,:]) - IRM_arr[i,i])*-1	
	for edge in tree.postorder_edge_iter():
		if edge.length != None:
			Pt_df = pd.DataFrame(scipy.linalg.expm(edge.length* IRM_arr), columns=all_states, index=all_states)
			Pt_dict[edge.head_node.taxon.label] = Pt_df
			

	
	for site in [str(x) for x in range(1, align_len+ 1)]:
		# remove gap character
		if options.exhaustive:
			obs_states = all_states
		elif options.constrained:
			# all rotamer states associated to observed AA states allows to consider rotamer states not observed at leaves for internal nodes
			list2d = [AA_2_rota(x) for x in constraint_split_align_df[site].unique() if x != '-']
			obs_states = list(itertools.chain(*list2d))
		else:
			obs_states = [x for x in split_align_df[site].unique() if x != '-']
		NODE_dict = {}

		# step 1 (leaves)
		for node in tree.leaf_node_iter():
			Cx_dict = {}
			Lx_dict = {}
			# current state is always the observed state
			current_state = split_align_df[site][node.taxon.label]
			# get the corresponding probability matrix
			Pt_df = Pt_dict[node.taxon.label]
			for parent_state in obs_states:
				if options.ignore_leaf == node.taxon.label:
					# treat leaf as internal node, allows to evaluate ancestral seq predictions on empirical alignments
					# where intnodes haven't been resurrected.
					# Done by ingnoring observed sequence, maximizing likelihood (Pt really as no childnodes) 
					# and assigning most likely current state for each parent state
					max_dict = {}
					for current_state in obs_states:
						# a leaf has no children nodes so Pt is the only factor 
						max_dict[current_state] = Pt_df[current_state][parent_state]
					max_current_state = max(max_dict, key=max_dict.get)
					maxL = max_dict[max_current_state]
					Cx_dict[parent_state] = max_current_state
					Lx_dict[parent_state] = maxL
				else:
					if current_state == '-':
						# set L to 1
						# this should carry no information upward
						Lx_dict[parent_state] = 1
					else:
						# TODO is this order arbitrary? (from is row, to is col)
						Lx_dict[parent_state] = Pt_df[current_state][parent_state]
					Cx_dict[parent_state] = current_state
			NODE_dict[node.taxon.label] = NODE(Cx_dict, Lx_dict)
		
		#each internal node visited after its children (postorder traversal)
		for node in tree.postorder_internal_node_iter():
			if node._parent_node is None:
				root = node
				# visit root, assign state (step 3)
				# root might have 2 or 3 children
				childNODEs = [NODE_dict[childnode.taxon.label] for childnode in root.child_node_iter()]
				max_dict = {}
				if options.constrained:
					constraint_state = constraint_split_align_df[site]['ROOT']
					for current_state in constrained_states:
						max_dict[current_state] = freqs_df[current_state].values[0] * np.prod(np.array([childNODE.Lx_dict[current_state] for childNODE in childNODEs]))
				else:
					# "collapsing AA reconstruction" approach: obs_states are collapsed roat states
					for current_state in obs_states:
						max_dict[current_state] = freqs_df[current_state].values[0] * np.prod(np.array([childNODE.Lx_dict[current_state] for childNODE in childNODEs]))
				#add a column to the reconstructed alignment	
				reconstructed_split_align_df[site] = ['NA']*len(reconstructed_split_align_df)
				reconstructed_split_align_df[site]['ROOT'] = max(max_dict, key=max_dict.get)
			
			else:
				# visit internal node (not root), perform step 2
				Cx_dict = {}
				Lx_dict = {}
				Pt_df = Pt_dict[node.taxon.label]
				# get the child NODEs
				childNODE_1, childNODE_2 = [NODE_dict[childnode.taxon.label] for childnode in node.child_node_iter()]
				
				for parent_state in obs_states:
					# max L | parent state
					max_dict = {}
					if options.constrained:
						# reconstructed AA state for internal node 
						# constraining rotamer state reconstruction
						constraint_state = constraint_split_align_df[site][node.taxon.label]
						# a subset of rotamer states corresponding to the constaint state
						constrained_states = AA_2_rota(constraint_state) 
						for current_state in constrained_states:
							max_dict[current_state] = Pt_df[current_state][parent_state] * childNODE_1.Lx_dict[current_state] * childNODE_2.Lx_dict[current_state]
					elif options.Ignore_and_trim_sibling in [childnode.taxon.label for childnode in node.child_node_iter()]:
						childNODE_1 = NODE_dict[options.Ignore_and_trim_sibling]
						for current_state in obs_states:
							# don't consider Lx from sibling of ignored node
							max_dict[current_state] = Pt_df[current_state][parent_state] * childNODE_1.Lx_dict[current_state]
					else:
						# for the "collapsing AA reconstruction" approach
						# here current state would actually be a tuple of rota states corresponding to one "current AA state"
						# so we are calculating the likelihood of I as the sum of I1, I2 amd I3
						# obs_states could be a list of collapsed observed states
						for current_state in obs_states:
							max_dict[current_state] = Pt_df[current_state][parent_state] * childNODE_1.Lx_dict[current_state] * childNODE_2.Lx_dict[current_state]
					max_current_state = max(max_dict, key=max_dict.get)
					maxL = max_dict[max_current_state]
					Cx_dict[parent_state] = max_current_state
					Lx_dict[parent_state] = maxL
				NODE_dict[node.taxon.label] = NODE(Cx_dict, Lx_dict)
			
		# assign states down the tree (step 5) to internal nodes based on state assigned to parent node
		nodes = [node for node in tree.preorder_internal_node_iter() if node is not root]
		if options.ignore_leaf:
			# also reconstruct ignored leaf
			nodes = nodes + [node for node in tree.leaf_node_iter() if options.ignore_leaf == node.taxon.label]
		for node in nodes:
			if node.parent_node is root:
				parent_assigned_state = reconstructed_split_align_df[site]['ROOT']
			else:
				parent_assigned_state = reconstructed_split_align_df[site][node.parent_node.taxon.label]
			current_assigned_state = NODE_dict[node.taxon.label].Cx_dict[parent_assigned_state]	
			reconstructed_split_align_df[site][node.taxon.label] = current_assigned_state
	
	reconstructed_split_align_df['sequences'] = reconstructed_split_align_df.iloc[:,0:align_len].apply(lambda x: ' '.join(x), axis=1)
	reconstructed_split_align_df['taxons'] = reconstructed_split_align_df.index

	print "%d\t%d" % (len(reconstructed_split_align_df), align_len)
	print
	if options.ignore_leaf:	
		# print only the reconstructed leaf sequence
		print "%s\t%s" % (options.ignore_leaf, reconstructed_split_align_df.loc[options.ignore_leaf]['sequences'])
	else:
		print reconstructed_split_align_df.to_csv(header=False, index=False, sep="\t", columns=['taxons', 'sequences'])

if __name__ == '__main__':
	main()

