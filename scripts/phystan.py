#!/usr/bin/env python
import argparse

import dendropy
import numpy
import pystan
import pickle
from dendropy import Tree, DnaCharacterMatrix
import os
import phylo

#
# reads a Stan ADVI diagnostic file to extract the best ELBO
#
def get_elbo(diag_filename):
    diag_file = open(diag_filename, "r")
    best_elbo=1
    for line in diag_file:
        if line.startswith("#"):
            continue
        line.rstrip("\n")
        d = line.split(",")
        d[2] = float(d[2])
        if best_elbo > 0 or best_elbo < d[2]:
            best_elbo = d[2]
    return best_elbo


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model', choices=['JC69', 'K80', 'HKY', 'GTR'], default='GTR',
                    help="""Substitution model [default: %(default)d]""")
parser.add_argument('-a', '--algorithm', choices=['vb', 'nuts', 'hmc'], default='vb', help="""Algoirthm [default: %(default)d]""")
parser.add_argument('-q', '--variational', choices=['meanfield', 'fullrank'], default='meanfield', help="""Variational approximation""")
parser.add_argument('-t', '--tree', type=argparse.FileType('r'), required=True, help="""Tree file""")
parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True, help="""Sequence file""")
parser.add_argument('-o', '--output', required=False, help="""Output file""")
parser.add_argument('-C', '--compile', action="store_true", help="""Recompile StanModel""")
arg = parser.parse_args()

my_path = os.path.split(os.path.realpath(__file__))[0]


taxa = dendropy.TaxonNamespace()

tree = Tree.get(file=arg.tree, schema="newick", tree_offset=0, taxon_namespace=taxa)

dna = DnaCharacterMatrix.get(file=arg.input, schema="fasta")
alignment_length = dna.sequence_size
sequence_count = len(dna)

print('Number of sequences: {} length {} '.format(sequence_count, alignment_length))
print('Model: ' + arg.model)

tipdata = phylo.get_dna_leaves_partials(dna)

for node in tree.postorder_node_iter():
    node.index = -1
    node.annotations.add_bound_attribute("index")

s = sequence_count + 1
for node in tree.postorder_node_iter():
    if not node.is_leaf():
        node.index = s
        s += 1
    else:
        for idx, name in enumerate(dna):
            if str(name) == str(node.taxon):
                node.index = idx + 1
                break

peeling = phylo.get_peeling_order(tree)

internal_nodes = []
nodes_height = {}
for node in tree.postorder_node_iter():
    if node.is_internal():
        child = node.child_node_iter().next()
        nodes_height[node] = nodes_height[child] + child.edge.length
        internal_nodes.append(node)
    else:
        nodes_height[node] = 0.0

internal_nodes.sort(key=lambda x: nodes_height[x], reverse=False)

# record lineage for each intercoalescent time (tau).
tau_map = numpy.zeros((sequence_count-1, sequence_count), dtype=numpy.int)

for idx, n in enumerate(internal_nodes):
    j = 0
    for node in tree.postorder_node_iter():
        if node.parent_node is not None and node != n and nodes_height[node.parent_node] >= nodes_height[n] >= nodes_height[node]:
            tau_map[idx][j] = node.index
            j += 1

data = {'peel': peeling, 'tipdata': tipdata, 'L': alignment_length, 'S': sequence_count, 'map': tau_map}

if arg.model == 'GTR':
    data['frequencies_alpha'] = [1, 1, 1, 1]
    data['rates_alpha'] = [1, 1, 1, 1, 1, 1]
    source = 'GTR-coalescent.stan'
elif arg.model == 'HKY':
    data['frequencies_alpha'] = [1, 1, 1, 1]
    source = 'HKY-coalescent.stan'
elif arg.model == 'K80':
    source = 'K80-coalescent.stan'
elif arg.model == 'JC69':
    source = 'JC69-coalescent.stan'

source_file = os.path.join(my_path, '..', 'src', source)

# Samples output file
sample_path = arg.input.name + '.log'
if arg.output:
    sample_path = arg.output

# Save binary to this file
binary = os.path.join(my_path, '..', 'bin', source + '.' + arg.algorithm + '.pkl')

# Save binary if it does not exists or compile option is set
if arg.compile or not os.path.isfile(binary):
    sm = pystan.StanModel(file=source_file)
    with open(binary, 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(binary, 'rb'))

if arg.algorithm == 'vb':
    fit = sm.vb(data=data, tol_rel_obj=0.001, elbo_samples=100, iter=10000, sample_file=sample_path, diagnostic_file=sample_path+".diag", algorithm=arg.variational)
    elbo = get_elbo(sample_path+".diag")
    print "Marginal likelihood lower bound is " + str(elbo)
else:
    fit = sm.sampling(data=data, algorithm=arg.algorithm.upper())
print(fit)

# with open(arg.input.name+'.info', 'w') as fp:
#     fp.write(str(fit))


