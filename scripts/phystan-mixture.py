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
parser.add_argument('-s', '--script', required=False, help="""Stan script file""")
parser.add_argument('-p', '--parameters', required=False, help="""Parameters for Stan script""")
parser.add_argument('-e', '--eta', required=False, type=int,  help="""Parameters for Stan script""")
arg = parser.parse_args()

my_path = os.path.split(os.path.realpath(__file__))[0]


taxa = dendropy.TaxonNamespace()

trees = dendropy.TreeList.get(file=arg.tree, schema="newick", preserve_underscores=True, tree_offset=0, taxon_namespace=taxa)

dna = DnaCharacterMatrix.get(file=arg.input, schema="fasta")
alignment_length = dna.sequence_size
sequence_count = len(dna)

print('Number of sequences: {} length {} '.format(sequence_count, alignment_length))
print('Model: ' + arg.model)


tipdata, weights = phylo.get_dna_leaves_partials_compressed(dna)
alignment_length = len(weights)

for t in trees:
    t.encode_bipartitions(collapse_unrooted_basal_bifurcation=False)

count = 1
bip = {}
bip2 = {}
bip3 = {}
indexes = []

# for tree in trees:
#     unwanted_node = None
#     for unwanted_node in tree.seed_node.child_node_iter():
#         pass
#     for node in tree.postorder_node_iter():
#         if node != tree.seed_node and unwanted_node != node:
#             if node.edge.bipartition not in bip2:
#                 bip2[node.edge.bipartition] = 1
#             else:
#                 bip2[node.edge.bipartition] += 1
#
# for tree in trees:
#     unwanted_node = None
#     for unwanted_node in tree.seed_node.child_node_iter():
#         pass
#     for node in tree.postorder_node_iter():
#         if node != tree.seed_node and unwanted_node != node:
#             if node.edge.bipartition not in bip2 or bip2[node.edge.bipartition] != len(trees):
#                 for x in node.adjacent_nodes():
#                     del bip2[x.edge.bipartition]
#
# for tree in trees:
#     unwanted_node = None
#     for unwanted_node in tree.seed_node.child_node_iter():
#         pass
#     for node in tree.postorder_node_iter():
#         if node.edge.bipartition in bip2:
#             count += 1
# print(count)
# exit(1)


# unwanted_node = None
# for unwanted_node in trees[0].seed_node.child_node_iter():
#     pass
#
# for node in trees[0].postorder_node_iter():
#     if node != trees[0].seed_node and unwanted_node != node:
#         bip[node.edge.bipartition] = count
#         indexes.append(count)
#         count += 1
# print(len(indexes))
#
#
# unwanted_node2 = None
# for unwanted_node2 in trees[1].seed_node.child_node_iter():
#     pass
#
# present = {}
# absent = {}
# for node in trees[1].postorder_node_iter():
#     if node != trees[1].seed_node and unwanted_node2 != node:
#         if node.edge.bipartition in bip:
#             present[node] = True
#         else:
#             present[node] = False
#
# for node in trees[1].postorder_node_iter():
#     if node != trees[1].seed_node and unwanted_node2 != node:
#         if not present[node]:
#             absent[node] = True
#             for x in node.adjacent_nodes():
#                 absent[x] = True
#         else:
#             absent[node] = False
#
# for node in trees[1].postorder_node_iter():
#     if node != trees[1].seed_node and unwanted_node2 != node:
#         if absent[node]:
#             indexes.append(count)
#             count += 1
#         else:
#             indexes.append(bip[node.edge.bipartition])
#
# print(indexes)
# print(max(indexes))


for tree in trees:
    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

leaf_indexes = {}
for idx, name in enumerate(dna):
    leaf_indexes[str(name)] = idx + 1

for tree in trees:
    s = sequence_count + 1
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            s += 1
        else:
            node.index = leaf_indexes[str(node.taxon)]

# peelings = phylo.get_peeling_orders(trees)
peelings = []
for tree in trees:
    peeling = phylo.get_peeling_order(tree)
    peelings.append(peeling)


coalescent = False
if coalescent:
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

if len(trees) == 1:
    data = {'peel': peelings[0], 'tipdata': tipdata, 'L': alignment_length, 'S': sequence_count}
else:
    data = {'peel': peelings, 'tipdata': tipdata, 'L': alignment_length, 'S': sequence_count, 'T': len(trees), 'weights': weights}
#print(peelings)


if len(trees) > 1:
    source = 'JC69-mixture.stan'
    #data['B'] = max(indexes)
    #data['indexes'] = indexes
elif arg.model == 'GTR':
    if coalescent:
        data['frequencies_alpha'] = [1, 1, 1, 1]
        data['rates_alpha'] = [1, 1, 1, 1, 1, 1]
        source = 'GTR-coalescent.stan'
        data['map'] = tau_map
    else:
        source = 'GTR-uniform.stan'
elif arg.model == 'HKY':
    data['frequencies_alpha'] = [1, 1, 1, 1]
    source = 'HKY-coalescent.stan'
    data['map'] = tau_map
elif arg.model == 'K80':
    source = 'K80-coalescent.stan'
    data['map'] = tau_map
elif arg.model == 'JC69':
    source = 'JC69-coalescent.stan'
    if coalescent:
        data['map'] = tau_map
    else:
        source = 'JC69-uniform.stan'

# script provided
if arg.script:
    source = arg.script
    source_file = source
    binary = source_file + '.pkl'
else:
    source_file = os.path.join(my_path, '..', 'src', source)
    # Save binary to this file
    binary = os.path.join(my_path, '..', 'bin', source + '.' + arg.algorithm + '.pkl')

# use arguments
# example: "frequencies_alpha=[1,1,1,1] rates_alpha=[1,1,1,1,1]"
if arg.parameters:
    params = arg.parameters.split(' ')
    for param in params:
        name, values = param.split('=')
        # list
        if values.startswith('['):
            data[name] = [float(v) for v in values.strip('[]').split()]
        # scalar
        else:
            data[name] = float(values)

# Samples output file
sample_path = arg.input.name + '.log'
if arg.output:
    sample_path = arg.output

# Save binary if it does not exists or compile option is set
if arg.compile or not os.path.isfile(binary):
    sm = pystan.StanModel(file=source_file)
    with open(binary, 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(binary, 'rb'))

if arg.algorithm == 'vb':
    kwargs = {}
    if arg.eta:
        kwargs['eta'] = arg.eta
        kwargs['adapt_engaged'] = False
    #fit = sm.vb(data=data, tol_rel_obj=0.000000001, elbo_samples=100, iter=1000000, sample_file=sample_path, diagnostic_file=sample_path+".diag", algorithm=arg.variational)
    fit = sm.vb(data=data, tol_rel_obj=0.001, elbo_samples=100, iter=10000, sample_file=sample_path,
                diagnostic_file=sample_path + ".diag", algorithm=arg.variational, **kwargs)
    #elbo = get_elbo(sample_path+".diag")
    #print "Marginal likelihood lower bound is " + str(elbo)
else:
    fit = sm.sampling(data=data, algorithm=arg.algorithm.upper())
#print(fit)

# with open(arg.input.name+'.info', 'w') as fp:
#     fp.write(str(fit))


