#!/usr/bin/env python
import argparse

import dendropy
import numpy
import pystan
import pickle
from dendropy import Tree, DnaCharacterMatrix
import os
import phylo
import generate_script
import glob
import tempfile
import filecmp
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


parser = argparse.ArgumentParser('Phylogenetics in Stan')
parser.add_argument('-m', '--model', choices=['JC69', 'GTR'], default='GTR',
                    help="""Substitution model [default: %(default)s]""")
parser.add_argument('-a', '--algorithm', choices=['vb', 'nuts', 'hmc'], default='vb',
                    help="""Algorithm [default: %(default)s]""")
parser.add_argument('-q', '--variational', choices=['meanfield', 'fullrank'], default='meanfield',
                    help="""Variational approximation""")
parser.add_argument('-t', '--tree', type=argparse.FileType('r'), required=True, help="""Tree file""")
parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True, help="""Sequence file""")
parser.add_argument('-o', '--output', required=False, help="""Output file""")
parser.add_argument('-F', '--force', action="store_true", help="""Recompile StanModel""")
parser.add_argument('-s', '--script', required=False, help="""Stan script file""")
parser.add_argument('-p', '--parameters', required=False, help="""Parameters for Stan script""")
parser.add_argument('-e', '--eta', required=False, type=float, help="""eta for Stan script""")
parser.add_argument('-S', '--seed', required=False, type=int, help="""Seed for Stan script""")
parser.add_argument('-I', '--invariant', required=False, action='store_true',
                    help="""Include a proportion of invariant sites""")
parser.add_argument('-C', '--categories', metavar='C', required=False, type=int, default=1,
                    help="""Number of categories""")
parser.add_argument('--heterogeneity', choices=['weibull', 'discrete'], default='weibull',
                    help="""Weibull or discrete distribution to model rate heterogeneity across sites""")
parser.add_argument('--heterochronous', action="store_true",
                    help="""Heterochronous data. Expect a date in the leaf names""")
parser.add_argument('--lower_root', type=float, default=0.0, help="""Lower bound of the root""")
parser.add_argument('--rate', required=False, type=float, help="""Susbstitution rate""")
parser.add_argument('--clock', required=False, choices=['strict'], default=None, help="""Type of clock""")
parser.add_argument('--estimate_rate', action='store_true', help="""Estimate substitution rate""")
parser.add_argument('-c', '--coalescent', choices=['constant', 'skyride', 'skygrid'], default=None,
                    help="""Type of coalescent (constant or skyride)""")
parser.add_argument('--grid', metavar='I', required=False, type=int, help="""Number of grid points in skygrid""")
parser.add_argument('--cutoff', metavar='G', required=False, type=float, help="""a cutoff for skygrid""")
arg = parser.parse_args()

my_path = os.path.split(os.path.realpath(__file__))[0]
bin_path = os.path.join(my_path, '..', 'bin')

if not os.path.lexists(bin_path):
    os.mkdir(bin_path)

taxa = dendropy.TaxonNamespace()

tree = Tree.get(file=arg.tree, schema="newick", tree_offset=0, taxon_namespace=taxa, preserve_underscores=True)

dna = DnaCharacterMatrix.get(file=arg.input, schema="nexus", preserve_underscores=True)
alignment_length = dna.sequence_size
sequence_count = len(dna)

print('Number of sequences: {} length {} '.format(sequence_count, alignment_length))
print('Model: ' + arg.model)


tipdata, weights = phylo.get_dna_leaves_partials_compressed(dna)
alignment_length = len(weights)

dates = []

for node in tree.postorder_node_iter():
    node.index = -1
    node.annotations.add_bound_attribute("index")

s = sequence_count + 1
for node in tree.postorder_node_iter():
    if not node.is_leaf():
        node.index = s
        s += 1
    else:
        found = False
        for idx, name in enumerate(dna):
            if str(name) == str(node.taxon):
                node.index = idx + 1
                found = True
                break
        if not found:
            print('Could not find taxon {} in alignment'.format(str(node.taxon)))
            exit(1)

# parse dates
if arg.heterochronous:
    newest = -1
    for node in tree.postorder_node_iter():
        node.date = 0.0
        node.annotations.add_bound_attribute("date")
        if node.is_leaf():
            node.date = float(str(node.taxon).split('_')[-1][:-1])
            newest = max(node.date, newest)
            dates.append(node.date)

    oldest = -1
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.date = newest - node.date
            oldest = max(oldest, node.date)

peeling = phylo.get_peeling_order(tree)

data = {'peel': peeling, 'tipdata': tipdata, 'L': alignment_length, 'S': sequence_count, 'weights': weights}

if arg.clock is not None:
    data['map'] = phylo.get_preorder(tree)
    data['rate'] = arg.rate if arg.rate else 1.0 # 7.9E-4
    if arg.heterochronous:
        data['lowers'] = phylo.get_lowers(tree)
        data['lower_root'] = max(oldest, arg.lower_root)
    else:
        data['lower_root'] = arg.lower_root

if arg.categories > 1:
    data['C'] = arg.categories
    if arg.invariant:
        data['C'] += 1

if arg.clock is not None:
    if arg.coalescent == 'skygrid':
        data['G'] = arg.grid - 1
        data['grid'] = numpy.linspace(0, arg.cutoff, arg.grid)[1:]
        # data['log_mu'] = numpy.log(4)

    # Not used by constant coalescent
    if arg.heterochronous:
        data['I'] = len(set(dates)) + sequence_count - 2 # #unique sampling dates + #internal nodes -1
    else:
        data['I'] = sequence_count - 1  # #unique sampling dates + #internal nodes -1

if arg.model == 'GTR':
    data['frequencies_alpha'] = [1, 1, 1, 1]
    data['rates_alpha'] = [1, 1, 1, 1, 1, 1]
elif arg.model == 'HKY':
    data['frequencies_alpha'] = [1, 1, 1, 1]

source_file = 'temp.stan'
binary = None
source_bin = None

# script provided
if arg.script:
    source_file = arg.script
    binary = source_file + '.pkl'
else:
    with open(source_file, 'w') as fp:
        script = generate_script.get_model(arg)
        fp.write(script)

stanscripts = glob.glob(os.path.join(my_path, '..', 'bin', '*.stan'))
for filename in stanscripts:
    f = os.path.join(my_path, '..', 'bin', filename)
    if filecmp.cmp(f, source_file, False):
        # stan file might exist but not the binary (script was killed)
        if os.path.isfile(f.replace('.stan', '.pkl')):
            binary = f.replace('.stan', '.pkl')
        source_bin = f
        break

# binary file does not exists so copy file to bin and then compile
if binary is None or arg.force:
    if binary is None and source_bin is None:
        f = tempfile.NamedTemporaryFile(dir=os.path.join(my_path, '..', 'bin'), suffix='.stan', delete=False)
        f.write(script)
        f.close()
        source_file = f.name
    else:
        source_file = source_bin
    print(source_file)
    sm = pystan.StanModel(file=source_file)
    binary = source_file.replace('.stan', '.pkl')
    with open(binary, 'wb') as f:
        pickle.dump(sm, f)
else:
    sm = pickle.load(open(binary, 'rb'))

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

if arg.algorithm == 'vb':
    stan_args = {}
    if arg.eta:
        stan_args['eta'] = arg.eta
        stan_args['adapt_engaged'] = False
    if arg.seed:
        stan_args['seed'] = arg.seed

    fit = sm.vb(data=data, tol_rel_obj=0.001, elbo_samples=100, iter=10000, sample_file=sample_path,
                diagnostic_file=sample_path + ".diag", algorithm=arg.variational, **stan_args)
    #elbo = get_elbo(sample_path+".diag")
    #print "Marginal likelihood lower bound is " + str(elbo)
else:
    fit = sm.sampling(data=data, algorithm=arg.algorithm.upper(), sample_file=sample_path)
#print(fit)

# with open(arg.input.name+'.info', 'w') as fp:
#     fp.write(str(fit))


