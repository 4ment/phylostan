#!/usr/bin/env python

import argparse
import os
import pickle
import sys

import dendropy
import numpy as np
import pystan
from dendropy import DnaCharacterMatrix, Tree

from . import utils
from .generate_script import get_model


def create_parse_parser(subprasers):
    parser = subprasers.add_parser('parse', help='parse Stan log files')
    parser.add_argument(
        '--samples', required=True, help="""Path to sample file from Stan"""
    )
    parser.add_argument(
        '-t', '--tree', type=argparse.FileType('r'), required=True, help="""Tree file"""
    )
    parser.add_argument('-o', '--output', required=True, help="""Nexus output file""")
    parser.add_argument(
        '--alpha',
        type=float,
        required=False,
        default=0.05,
        help="""Controls level for 100*(1-alpha)% Bayesian credible intervals""",
    )
    parser.add_argument(
        '--rate', type=float, required=False, help="""Value of fixed rate"""
    )
    parser.add_argument(
        '--dates',
        required=False,
        type=argparse.FileType('r'),
        help="""Comma-separated (csv) file containing sequence dates with
        header 'name,date'""",
    )
    parser.add_argument(
        '--heterochronous',
        action="store_true",
        help="""Heterochronous data. Expect a date in the leaf names or a csv file
        containing dates""",
    )
    return parser


def create_run_parser(subprasers):
    parser = create_build_parser(
        subprasers, 'run', help='run an analysis using a Stan script'
    )
    parser.add_argument(
        '-t', '--tree', type=argparse.FileType('r'), required=True, help="""Tree file"""
    )
    parser.add_argument(
        '-i',
        '--input',
        type=argparse.FileType('r'),
        required=False,
        help="""Sequence file""",
    )
    parser.add_argument(
        '-o', '--output', required=True, help="""Stem for output files"""
    )
    parser.add_argument(
        '--lower_root', type=float, default=0.0, help="""Lower bound of the root"""
    )
    parser.add_argument(
        '--rate', required=False, type=float, help="""Substitution rate"""
    )
    parser.add_argument(
        '--dates',
        required=False,
        type=argparse.FileType('r'),
        help="""Comma-separated (csv) file containing sequence dates with header
        'name,date'""",
    )
    parser.add_argument(
        '--heights_init',
        action="store_true",
        help="""initialize node heights using input tree file""",
    )

    parser.add_argument(
        '-a',
        '--algorithm',
        choices=['VB', 'NUTS', 'HMC', 'LBFGS'],
        default='VB',
        type=str.upper,
        help="""Algorithm [default: %(default)s]""",
    )
    parser.add_argument(
        '-S', '--seed', required=False, type=int, help="""Seed for Stan script"""
    )
    # variational
    parser.add_argument(
        '-q',
        '--variational',
        choices=['meanfield', 'fullrank'],
        default='meanfield',
        help="""Variational distribution family""",
    )
    parser.add_argument(
        '-e',
        '--eta',
        required=False,
        type=float,
        help="""eta for Stan script (variational only)""",
    )
    parser.add_argument(
        '--elbo_samples',
        required=False,
        type=int,
        default=100,
        help="""Number of samples for Monte Carlo estimate of ELBO
        (variational only)""",
    )
    parser.add_argument(
        '--grad_samples',
        required=False,
        type=int,
        default=1,
        help="""Number of samples for Monte Carlo estimate of gradients
    (variational only)""",
    )
    parser.add_argument(
        '--samples',
        required=False,
        type=int,
        default=1000,
        help="""Number of samples to be drawn from the variational distribution
        (variational only)""",
    )
    parser.add_argument(
        '--tol_rel_obj',
        required=False,
        type=float,
        default=0.001,
        help="""Convergence tolerance on the relative norm of the objective, defaults
        to 0.001 (variational only)""",
    )
    # sampling
    parser.add_argument(
        '--chains',
        required=False,
        type=int,
        default=1,
        help="""Number of chains for sampling algorithms (NUTS and HMC)""",
    )
    parser.add_argument(
        '--thin',
        required=False,
        type=int,
        default=1,
        help="""Positive integer specifying the period for saving samples
        (NUTS and HMC)""",
    )
    # sampling and variational
    parser.add_argument(
        '--iter',
        required=False,
        type=int,
        default=100000,
        help="""Maximum number of iterations for variational inference or number
        of iterations for NUTS and HMC algorithms""",
    )
    parser.add_argument('-M', '--metadata', help="""Phylogeography metadata file""")
    parser.add_argument(
        '--metadata_key',
        default="Country",
        help="""Phylogeography: use key in metadata""",
    )

    parser.add_argument(
        '--init',
        type=argparse.FileType('r'),
        help="""File containing starting values""",
    )
    return parser


def create_build_parser(subprasers, prog, help):
    parser = subprasers.add_parser(prog, help=help)
    parser.add_argument('-s', '--script', required=True, help="""Stan script file""")
    parser.add_argument(
        '-m',
        '--model',
        choices=['JC69', 'HKY', 'GTR'],
        default='GTR',
        help="""Substitution model [default: %(default)s]""",
    )
    parser.add_argument(
        '-I',
        '--invariant',
        required=False,
        action='store_true',
        help="""Include a proportion of invariant sites""",
    )
    parser.add_argument(
        '-C',
        '--categories',
        metavar='C',
        required=False,
        type=int,
        default=1,
        help="""Number of categories""",
    )
    parser.add_argument(
        '--heterogeneity',
        choices=['weibull', 'discrete'],
        default='weibull',
        help="""Weibull or discrete distribution to model rate heterogeneity
        across sites""",
    )
    parser.add_argument(
        '--heterochronous',
        action="store_true",
        help="""Heterochronous data. Expect a date in the leaf names""",
    )
    parser.add_argument(
        '--clock',
        required=False,
        choices=[
            'strict',
            'ace',
            'acln',
            'acg',
            'aoup',
            'ucln',
            'uced',
            'gmrf',
            'hsmrf',
        ],
        default=None,
        help="""Type of clock""",
    )
    parser.add_argument(
        '--clockpr',
        default='ctmcscale',
        type=lambda x: distribution_type(x, ('exponential', 'ctmcscale')),
        help="""prior on substitution rate [default: %(default)s]""",
    )
    parser.add_argument(
        '--estimate_rate', action='store_true', help="""Estimate substitution rate"""
    )
    parser.add_argument(
        '-c',
        '--coalescent',
        choices=['constant', 'skyride', 'skygrid', 'skyglide'],
        default=None,
        help="""Type of coalescent (constant or skyride)""",
    )
    parser.add_argument(
        '--speciation',
        choices=['bd', 'yule'],
        default=None,
        help="""Speciation model (birth-death or Yule)""",
    )
    parser.add_argument(
        '--grid',
        metavar='I',
        required=False,
        type=int,
        help="""Number of grid points in skygrid or skyglide""",
    )
    parser.add_argument(
        '--cutoff',
        metavar='G',
        required=False,
        type=float,
        help="""a cutoff for skygrid and skyglide""",
    )
    parser.add_argument(
        '--non_centered',
        action="store_true",
        help="""Use non centered parameterization for population size parameters
        (skygrid and skyglide only)""",
    )
    parser.add_argument(
        '--time_aware',
        action="store_true",
        help="""Use time-aware GMRF (skyride only)""",
    )
    parser.add_argument(
        '--compile', action="store_true", help="""Compile Stan script"""
    )
    parser.add_argument(
        '--rescaling',
        action="store_true",
        help="""Use rescaling in DNA tree likelihood calculation""",
    )
    parser.add_argument('--geo', action="store_true", help="""Phylogeography""")
    parser.add_argument(
        '--rescaling_geo', action="store_true", help="""Use rescaling phylogeography"""
    )
    return parser


def create_compile_parser(subprasers, prog, help):
    parser = subprasers.add_parser(prog, help=help)
    parser.add_argument('-s', '--script', required=True, help="""Stan script file""")
    return parser


def distribution_type(arg, choices):
    """Used by argparse for specifying distributions with optional
    parameters."""
    res = arg.split('(')
    if (isinstance(choices, tuple) and res[0] in choices) or res[0] == choices:
        return arg
    else:
        if isinstance(choices, tuple):
            message = "'" + "','".join(choices) + '"'
        else:
            message = "'" + choices + "'"
        raise argparse.ArgumentTypeError(
            'invalid choice (choose from a number or ' + message + ')'
        )


def main():
    parser_epilog = """To get some information for each sub-command:\n
    phylostan build --help
    phylostan run --help
    phylostan parse --help
    phylostan compile --help
"""

    parser = argparse.ArgumentParser(
        prog='phylostan',
        description='Phylogenetic inference using Stan',
        epilog=parser_epilog,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subprasers = parser.add_subparsers()

    build_parser = create_build_parser(subprasers, 'build', 'build a Stan script')
    build_parser.set_defaults(func=build)

    run_parser = create_run_parser(subprasers)
    run_parser.set_defaults(func=run)

    parse_parser = create_parse_parser(subprasers)
    parse_parser.set_defaults(func=parse)

    compile_parser = create_compile_parser(
        subprasers, 'compile', 'compile a Stan script'
    )
    compile_parser.set_defaults(func=compile_script)

    # Set default function to show help if no subcommand is used
    def show_help(arg):
        parser.print_help()
        sys.exit(1)

    parser.set_defaults(func=show_help)

    arg = parser.parse_args()
    arg.func(arg)


def parse_logs(treeobj, treelog, samplelog, rate, alpha):
    utils.convert_samples_to_nexus(treeobj, samplelog, treelog, rate)
    utils.parse_log(samplelog, alpha, treeobj)


def parse(arg):
    taxa = dendropy.TaxonNamespace()

    tree_format = 'newick'
    with open(arg.tree.name) as fp:
        if next(fp).upper().startswith('#NEXUS'):
            tree_format = 'nexus'

    tree = Tree.get(
        file=arg.tree,
        schema=tree_format,
        tree_offset=0,
        taxon_namespace=taxa,
        preserve_underscores=True,
        rooting='force-rooted',
    )
    if len(tree.seed_node.adjacent_nodes()) > 2:
        tree.reroot_at_edge(tree.seed_node.adjacent_nodes()[0].edge)

    utils.setup_indexes(tree)
    utils.setup_dates(tree, arg.dates, arg.heterochronous)

    parse_logs(tree, arg.output, arg.samples, arg.rate, arg.alpha)


def build(arg):
    with open(arg.script, 'w') as fp:
        script = get_model(arg)
        fp.write(script)

    if arg.compile:
        compile_script_func(arg.script)


def compile_script(arg):
    compile_script_func(arg.script)


def compile_script_func(script_name):
    binary = script_name.replace('.stan', '.pkl')
    # arg.script does not end with .stan
    if binary == script_name:
        binary = script_name + '.pkl'
    sm = pystan.StanModel(file=script_name)
    with open(binary, 'wb') as f:
        pickle.dump(sm, f)


def run(arg):
    taxa = dendropy.TaxonNamespace()

    tree_format = 'newick'
    with open(arg.tree.name) as fp:
        if next(fp).upper().startswith('#NEXUS'):
            tree_format = 'nexus'

    tree = Tree.get(
        file=arg.tree,
        schema=tree_format,
        tree_offset=0,
        taxon_namespace=taxa,
        preserve_underscores=True,
        rooting='force-rooted',
    )

    tree.resolve_polytomies(update_bipartitions=True)

    utils.setup_indexes(tree)

    oldest = utils.setup_dates(tree, arg.dates, arg.heterochronous)

    peeling = utils.get_peeling_order(tree)
    sequence_count = len(tree.taxon_namespace)
    data = {'peel': peeling, 'S': sequence_count}

    if arg.input:
        seqs_args = dict(schema='nexus', preserve_underscores=True)
        with open(arg.input.name) as fp:
            if next(fp).startswith('>'):
                seqs_args = dict(schema='fasta')

        dna = DnaCharacterMatrix.get(file=arg.input, taxon_namespace=taxa, **seqs_args)
        alignment_length = dna.sequence_size
        sequence_count = len(dna)
        if sequence_count != len(dna.taxon_namespace):
            sys.stderr.write('taxon names in trees and alignment are different')
            exit(2)

        print(
            'Number of sequences: {} length {} '.format(
                sequence_count, alignment_length
            )
        )
        print('Model: ' + arg.model)

        tipdata, weights = utils.get_dna_leaves_partials_compressed(dna)
        alignment_length = len(weights)

        data.update({'tipdata': tipdata, 'L': alignment_length, 'weights': weights})

    if arg.metadata:
        # Parse metadata file
        with open(arg.metadata) as fp:
            geodata = {}
            countries = {}
            geopattern = []
            header = next(fp).strip().split('\t')
            index_country = header.index(arg.metadata_key)

            for line in fp:
                row = line.strip().split('\t')
                if len(row) > 0:
                    geodata[row[0]] = row[index_country]
                    countries[row[index_country]] = 1

        country_to_index = {}
        index_to_country = []
        for idx, taxon in enumerate(tree.taxon_namespace):
            country = geodata[taxon.label]
            if country not in country_to_index:
                country_to_index[country] = len(country_to_index)
                index_to_country.append(country)

        print('"' + '","'.join(index_to_country) + '"')

        state_count = len(country_to_index)

        for idx, taxon in enumerate(tree.taxon_namespace):
            pattern = [0] * state_count
            country = geodata[taxon.label]
            pattern[country_to_index[country]] = 1
            geopattern.append(pattern)

        blens = [None] * (sequence_count * 2 - 1)
        for node in tree.postorder_node_iter():
            blens[node.index - 1] = node.edge.length
            if node.edge.length < 0:
                exit(3)

        children = tree.seed_node.child_nodes()
        blens[children[0].index] += blens[children[1].index]
        blens = blens[:-2]  # discard root branch and one of its child

        data['STATES'] = state_count
        data['blens'] = blens
        data['frequencies_alpha_geo'] = [1] * state_count
        data['rates_alpha_geo'] = [1] * int(state_count * (state_count - 1) / 2)
        data['geodata'] = geopattern

    if arg.clock is not None:
        data['map'] = utils.get_preorder(tree)
        if not arg.estimate_rate:
            data['rate'] = arg.rate if arg.rate else 1.0
        if arg.heterochronous:
            data['lowers'] = utils.get_lowers(tree)
            data['lower_root'] = max(oldest, arg.lower_root)
        else:
            data['lower_root'] = arg.lower_root
    else:
        last = peeling[-1]
        if last[0] > last[1]:
            peeling[-1] = [last[1], last[0], last[2]]

    if arg.categories > 1:
        data['C'] = arg.categories
        if arg.invariant:
            data['C'] += 1

    if arg.clock is not None:
        if arg.coalescent in ('skygrid', 'skyglide'):
            data['G'] = arg.grid - 1
            data['grid'] = np.linspace(0, arg.cutoff, arg.grid)[1:]
        elif arg.coalescent == 'skyride':
            # number of coalescent intervals
            data['I'] = sequence_count - 1

    if arg.model == 'GTR':
        data['frequencies_alpha'] = [1, 1, 1, 1]
        data['rates_alpha'] = [1, 1, 1, 1, 1, 1]
    elif arg.model == 'HKY':
        data['frequencies_alpha'] = [1, 1, 1, 1]

        # Samples output file
    sample_path = arg.output
    tree_path = sample_path + '.trees'

    binary = arg.script.replace('.stan', '.pkl')
    if binary == arg.script:
        binary = arg.script + '.pkl'
    if not os.path.lexists(binary) or arg.compile:
        sm = pystan.StanModel(file=arg.script)
        with open(binary, 'wb') as f:
            pickle.dump(sm, f)
    else:
        sm = pickle.load(open(binary, 'rb'))

    stan_args = {
        'data': data,
        'iter': arg.iter,
        'sample_file': sample_path,
        'algorithm': arg.algorithm,
    }
    if hasattr(arg, 'seed'):
        stan_args['seed'] = arg.seed

    if arg.init is not None:
        inits = {}
        for line in arg.init:
            line = line.strip()
            if len(row) > 0:
                line = line.split(':')
                inits[line[0].strip()] = list(map(float, line[1].split(',')))
        stan_args['init'] = inits
    elif arg.heights_init or arg.rate is not None:
        inits = {}
        if arg.heights_init:
            ratios, root_height = utils.ratios_root_height_from_branch_lengths(tree)
            # ratios_unres = np.log(ratios / (1.0 - ratios))
            # root_height_unres = np.log(root_height - data['lower_root'])
            inits['props'] = ratios.tolist()  # ratios_unres.tolist()
            inits['height'] = root_height.item() - data['lower_root']
            inits['rate'] = arg.rate
        elif arg.rate is not None:
            inits['rate'] = arg.rate
        stan_args['init'] = inits

    if arg.algorithm == 'LBFGS':
        fit = sm.optimizing(**stan_args)
        print(fit)
    elif arg.algorithm == 'VB':
        stan_args['algorithm'] = arg.variational
        stan_args['output_samples'] = arg.samples
        if arg.eta:
            stan_args['eta'] = arg.eta
            stan_args['adapt_engaged'] = False

        fit = sm.vb(
            tol_rel_obj=arg.tol_rel_obj,
            elbo_samples=arg.elbo_samples,
            grad_samples=arg.grad_samples,
            diagnostic_file=sample_path + ".diag",
            **stan_args
        )

        # parse the log file
        utils.convert_samples_to_nexus(tree, sample_path, tree_path, arg.rate)
        utils.parse_log(sample_path, 0.05)
    else:
        fit = sm.sampling(chains=arg.chains, thin=arg.thin, **stan_args)

        # chain=1 pystan uses sample_file
        if arg.chains == 1:
            if sample_path.endswith('.csv'):
                tree_path = sample_path.replace('.csv', '.trees')
            utils.convert_samples_to_nexus(tree, sample_path, tree_path, arg.rate)
            utils.parse_log(sample_path, 0.05)
            # chain>1 pystan appends _{chain}.csv to sample_file
        else:
            for chain in range(arg.chains):
                sample_path_chain = sample_path + '_{}.csv'.format(chain)
                tree_path_chain = sample_path + '_{}.trees'.format(chain)
                utils.convert_samples_to_nexus(
                    tree, sample_path_chain, tree_path_chain, arg.rate
                )
                utils.parse_log(sample_path_chain, 0.05)


if __name__ == "__main__":
    main()
