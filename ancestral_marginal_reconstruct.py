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


one_rota_to_four_dict = {'1': 'VAL2', '0': 'VAL1', '2': 'VAL3', 'A': 'ALA', 'C': 'ARG2', 'B': 'ARG1', 'E': 'ASN1', 'D': 'ARG3', 'G': 'ASN3', 'F': 'ASN2', 'I': 'ASP2', 'H': 'ASP1', 'K': 'CYS1', 'J': 'ASP3', 'M': 'CYS3', 'L': 'CYS2', 'O': 'GLN2', 'N': 'GLN1', 'Q': 'GLU1', 'P': 'GLN3', 'S': 'GLU3', 'R': 'GLU2', 'U': 'HIS1', 'T': 'GLY', 'W': 'HIS3', 'V': 'HIS2', 'Y': 'ILE2', 'X': 'ILE1', 'Z': 'ILE3', 'a': 'LEU1', 'c': 'LEU3', 'b': 'LEU2', 'e': 'LYS2', 'd': 'LYS1', 'g': 'MET1', 'f': 'LYS3', 'i': 'MET3', 'h': 'MET2', 'k': 'PHE2', 'j': 'PHE1', 'm': 'PRO1', 'l': 'PHE3', 'o': 'SER1', 'n': 'PRO2', 'q': 'SER3', 'p': 'SER2', 's': 'THR2', 'r': 'THR1', 'u': 'TRP1', 't': 'THR3', 'w': 'TRP3', 'v': 'TRP2', 'y': 'TYR2', 'x': 'TYR1', 'z': 'TYR3'}

four_rota_to_one_dict = dict([[v,k] for k,v in one_rota_to_four_dict.items()])

three_AA_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'TRP': 'W', 'PRO': 'P', 'THR': 'T', 'ILE': 'I', 'ALA': 'A', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'MET': 'M', 'GLU': 'E', 'ASN': 'N', 'TYR': 'Y', 'VAL': 'V'}

one_AA_to_three_dict = dict([[v,k] for k,v in three_AA_to_one.items()])

three_AA_to_three_rota_dict = {'CYS': ['CYS2', 'CYS1', 'CYS3'], 'ASP': ['ASP1', 'ASP3', 'ASP2'], 'SER': ['SER2', 'SER1', 'SER3'], 'GLN': ['GLN1', 'GLN3', 'GLN2'], 'LYS': ['LYS1', 'LYS3', 'LYS2'], 'ILE': ['ILE1', 'ILE3', 'ILE2'], 'PRO': ['PRO2', 'PRO1'], 'THR': ['THR1', 'THR3', 'THR2'], 'PHE': ['PHE1', 'PHE3', 'PHE2'], 'ALA': ['ALA'], 'GLY': ['GLY'], 'HIS': ['HIS2', 'HIS1', 'HIS3'], 'GLU': ['GLU2', 'GLU1', 'GLU3'], 'LEU': ['LEU2', 'LEU1', 'LEU3'], 'ARG': ['ARG1', 'ARG3', 'ARG2'], 'TRP': ['TRP2', 'TRP1', 'TRP3'], 'VAL': ['VAL1', 'VAL3', 'VAL2'], 'ASN': ['ASN2', 'ASN1', 'ASN3'], 'TYR': ['TYR1', 'TYR3', 'TYR2'], 'MET': ['MET2', 'MET1', 'MET3']}

four_rota_to_three_AA_dict = {}
for key, value in three_AA_to_three_rota_dict.items():
	for s in value:
		four_rota_to_three_AA_dict.setdefault(s, []).append(key)

def AA_2_rota(one_char_AA_state):
	rota_states = three_AA_to_three_rota_dict[one_AA_to_three_dict[one_char_AA_state]]
	return [four_rota_to_one_dict[r] for r in rota_states]

def rota_2_AA(one_char_rota_state):
	AA_states = four_rota_to_three_AA_dict[one_rota_to_four_dict[one_char_rota_state]]
	return [three_AA_to_one[a] for a in AA_states]

def main():
	usage = format_usage('''
		%prog NEWICK_TREE PHYLIP_ALIGN(rotamer states) IRM_SCALED FREQ_FILE

		
		''')
	parser = OptionParser(usage=usage)
	
	parser.add_option('-s', '--state_set', type=str, help='specify state set if IRM or freqs labels  and alignment symbols dont match')
	parser.add_option('-e', '--exhaustive', action="store_true", default=False, help='Iter over all states not only those observed at site')
	parser.add_option('-c', '--constrained', default=False, help='Constrain rotamer states reconstruction using internal AA sequences.')
	parser.add_option('-I', '--Ignore_and_trim_sibling', type=str, default=False, help='Ignore leaf sequences (expects file with ,-separated leafnames)')
	parser.add_option('-o', '--outfile', type=str, default=False, help='print to file as workaround for bsub')

	options, args = parser.parse_args()
	
	if len(args) != 4:
		exit('Unexpected argument number.')
	
	treefile, alignfile, IRMfile, freqfile = args

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
	
	def compute_assignments(reconstructed_split_align_df, ignore_leaf=False, trim_sibling=[]):
		internal_taxons = reconstructed_split_align_df.index.tolist()
		for site in [str(x) for x in range(1, align_len+ 1)]:
			P_dict = {}
			assign_dicts = []

			# the observed (1-character) ROTAMER STATES at the tips for this site
			obs_states = [x for x in split_align_df[site].unique() if x != '-']
			
			for assigned_states in itertools.combinations_with_replacement(obs_states, len(internal_taxons)):
				NODE_dict = dict(zip(internal_taxons, assigned_states))
				P_list = []
				# leaves
				for node in tree.leaf_node_iter():
					if node.taxon.label == ignore_leaf:
						# treat as internal node
						current_assigned_state = NODE_dict[node.taxon.label]
						Pt_df = Pt_dict[node.taxon.label]
						if node.parent_node.taxon != None:
							parent_assigned_state = NODE_dict[node.parent_node.taxon.label]
						else:
							parent_assigned_state = NODE_dict['ROOT']
						P_list.append(Pt_df[parent_assigned_state][current_assigned_state])
					elif node.taxon.label in trim_sibling:
						# ignore contribution from this node
						continue
					else:
						if split_align_df[site][node.taxon.label] == '-':
							# skip this taxon
							continue
						else:
							# current ROTAMER state is always the observed state
							current_state = split_align_df[site][node.taxon.label]
						# get the corresponding probability matrix
						Pt_df = Pt_dict[node.taxon.label]
					
						# get the parent node's assigned state
						parent = node.parent_node
						if parent.parent_node is None:
							parent_assigned_state = NODE_dict['ROOT']
							P_list.append(Pt_df[parent_assigned_state][current_state] * freqs_df[parent_assigned_state])
						else:
							parent_assigned_state = NODE_dict[parent.taxon.label]
					
							# get the transition probability from assigned parent state to observed tip state
							P_list.append(Pt_df[parent_assigned_state][current_state])
					
				# internal nodes including ROOT
				for node in tree.postorder_internal_node_iter():
					if node._parent_node is None:
						root = node
						current_assigned_state = NODE_dict['ROOT']
						for childnode in root.child_node_iter():
							# not dry
							if childnode.is_internal():
								Pt_df = Pt_dict[childnode.taxon.label]
								childnode_assigned_state = NODE_dict[childnode.taxon.label]
								P_list.append(Pt_df[current_assigned_state][childnode_assigned_state] * freqs_df[current_assigned_state])

					
					else:
						current_assigned_state = NODE_dict[node.taxon.label]
						for childnode in node.child_node_iter():
							if childnode.is_internal():
								Pt_df = Pt_dict[childnode.taxon.label]
								childnode_assigned_state = NODE_dict[childnode.taxon.label]
								P_list.append(Pt_df[current_assigned_state][childnode_assigned_state])
				# add the probability of this assignment to NODE_dict
				NODE_dict['P_assignment'] = (float(np.prod(P_list)))
				assign_dicts.append(NODE_dict)


			# add a column to the reconstructed alignment
			reconstructed_split_align_df[site] = ['NA']*len(reconstructed_split_align_df)

			# performs marginal reconstruction:
			# assign the MOST LIKELY STATE to each node INDEPENDETLY of other nodes' assignments
			for nodename in internal_taxons:
				best_state_dict = {}
				for obs_state in obs_states:
					# sum mrginal P for all assignmets that have obs_state at node
					s = sum([NODE_dict['P_assignment'] for NODE_dict in assign_dicts if NODE_dict[nodename] == obs_state])
					best_state_dict[obs_state] = s
				# get the assignment that corresponds to  the highest marginal P for this node (fastest method)
				v=list(best_state_dict.values())
				k=list(best_state_dict.keys())
				reconstructed_split_align_df[site][nodename] = k[v.index(max(v))]
		reconstructed_split_align_df['sequences'] = reconstructed_split_align_df.iloc[:,0:align_len].apply(lambda x: ' '.join(x), axis=1)
		reconstructed_split_align_df['taxons'] = reconstructed_split_align_df.index
		return reconstructed_split_align_df



	f_out = open(options.outfile, "a")
	if options.Ignore_and_trim_sibling:
		# iterate over the leaves in file
		for leafname in  open(options.Ignore_and_trim_sibling, "r").readlines()[0].split(","): 	
			in_df = pd.DataFrame(0, index=internal_taxons + [leafname], columns=[])
			ignore_leaf = leafname
			ignored_node = tree.find_node_with_taxon_label(leafname)

			# find the ignored mode's sibling
			trim_sibling = [node.taxon.label for node in ignored_node.sibling_nodes()]
		
			out_df = compute_assignments(in_df, ignore_leaf=ignore_leaf, trim_sibling=trim_sibling)	

			# print only the reconstructed leaf sequence
			print >>f_out, "%s\t%s" % (ignore_leaf, out_df.loc[ignore_leaf]['sequences'])
	else:
		in_df = pd.DataFrame(0, index=internal_taxons, columns=[])
		out_df = compute_assignments(in_df)
		print >>f_out, "%d\t%d" % (len(out_df), align_len)
		print >>f_out, ""
		print >>f_out, out_df.to_csv(header=False, index=False, sep="\t", columns=['taxons', 'sequences'])

if __name__ == '__main__':
	main()

