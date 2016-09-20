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
parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True, help="""Sequence file""")
parser.add_argument('-o', '--output', required=False, help="""Output file""")
parser.add_argument('-C', '--compile', action="store_true", help="""Recompile StanModel""")
arg = parser.parse_args()

my_path = os.path.split(os.path.realpath(__file__))[0]



dna = DnaCharacterMatrix.get(file=arg.input, schema="fasta")
alignment_length = dna.sequence_size
sequence_count = len(dna)


print('Number of sequences: {} length {} '.format(sequence_count, alignment_length))
print('Model: ' + arg.model)


##
# Prepare the Stan program
#

if arg.model == 'GTR':
    source = 'GTR-coalescent.stan'
elif arg.model == 'HKY':
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
binary = os.path.join(my_path, '..', 'bin', source + '.vb.pkl')

# Save binary if it does not exist or compile option is set
if arg.compile or not os.path.isfile(binary):
    sm = pystan.StanModel(file=source_file)
    with open(binary, 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(binary, 'rb'))

##
# begin the SMC.
# Add sequences in the order they appear in the input file
# Choose branches for addition uniformly at random
#
tipdata = phylo.get_dna_leaves_partials(dna)


# make the initial particles with all taxa from the first three sequences
n0 = dendropy.datamodel.treemodel.Node(label="0",edge_length=1)
n1 = dendropy.datamodel.treemodel.Node(label="1",edge_length=1)

t0n2 = dendropy.datamodel.treemodel.Node(label="n0",edge_length=1)
t0n2.add_child(n0)
t0n2.add_child(n1)
t0 = Tree(seed_node=t0n2)

# stores the current set of tree topology particles
psize=1000
particles = [t0]*psize

for i in range(2,sequence_count):
    # propose the addition of sequence i to particle p
    for p in particles:
        ni = dendropy.datamodel.treemodel.Node(label=str(i),edge_length=1)
        # pick an node at random
        nlist = p.nodes()
        n = nlist[numpy.random.randint(len(nlist))]
        # make a new internal node with n and ni as its child
        newnode = dendropy.datamodel.treemodel.Node(label="n" + str(i-1),edge_length=1)
        newnode.add_child(ni)
        parent = n.parent_node
        if parent is not None:
            parent.remove_child(n)
            newnode.add_child(n)
            parent.add_child(newnode)
        else:
            newnode.add_child(n)
            p.reroot_at_node(newnode)

    # calculate ELBOs for each unique topology
    cur_trees = set()
    for p in particles:
        cur_trees.add(p)
    particle_weights = dict()

    for tree in cur_trees:
        for node in tree.postorder_node_iter():
            node.index = -1
            node.annotations.add_bound_attribute("index")

        s = i + 1
        for node in tree.postorder_node_iter():
            if not node.is_leaf():
                node.index = s
                s += 1
            else:
                for idx, name in enumerate(dna):
                    if str(idx) == str(node.taxon):
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
        tau_map = numpy.zeros((i, i+1), dtype=numpy.int)

        for idx, n in enumerate(internal_nodes):
            j = 0
            for node in tree.postorder_node_iter():
                if node.parent_node is not None and node != n and nodes_height[node.parent_node] >= nodes_height[n] >= nodes_height[node]:
                    tau_map[idx][j] = node.index
                    j += 1

        data = {'peel': peeling, 'tipdata': tipdata, 'L': alignment_length, 'S': i+1, 'map': tau_map}


        if arg.model == 'GTR':
            data['frequencies_alpha'] = [1, 1, 1, 1]
            data['rates_alpha'] = [1, 1, 1, 1, 1, 1]
        elif arg.model == 'HKY':
            data['frequencies_alpha'] = [1, 1, 1, 1]


        fit = sm.vb(data=data, tol_rel_obj=0.001, elbo_samples=100, iter=10000, sample_file=sample_path, diagnostic_file=sample_path+".diag", algorithm=arg.variational)
        elbo = get_elbo(sample_path+".diag")

        # store the ELBO for this tree
        particle_weights[tree] = elbo

    # now resample the particles based on weight

