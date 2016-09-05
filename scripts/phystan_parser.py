#!/usr/bin/python

import sys
import dendropy
import phylo
import scipy
import argparse


def upadate_tree(tree, params):
    """Update branch lengths of tree with intercoalescent times parameters in params"""

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

    #taus.1,taus.2,taus.3,taus.4,taus.5,freqs.1,freqs.2,freqs.3,freqs.4,rates.1,rates.2,rates.3,rates.4,rates.5,rates.6,theta

    nodes_length = {}
    for idx, n in enumerate(internal_nodes):
        for node in tree.postorder_node_iter():
            if node.parent_node is not None and n != node and nodes_height[node.parent_node] >= nodes_height[n] >= nodes_height[node]:
                try:
                    nodes_length[node] += float(params['taus.{}'.format(idx+1)])
                except KeyError:
                    nodes_length[node] = float(params['taus.{}'.format(idx+1)])

    for node in tree.postorder_node_iter():
        if node.parent_node is not None:
            node.edge.length = nodes_length[node]
        else:
            node.edge.length = 0.0
    return tree


def recompute_all(seq_path, trees_path):
    dna = dendropy.DnaCharacterMatrix.get(path=seq_path, schema="fasta")
    taxa = dendropy.TaxonNamespace()

    output_probs = open(seq_path + '.csv', "w")
    output_probs.write('sample\ttree\tss\tposterior\tlikelihood\tprior\tCoalescentConstant\trootHeight\tfreqParameter.1\tfreqParameter.2\tfreqParameter.3\tfreqParameter.4')
    output_probs.write('\trateAC\trateAG\trateAT\trateCG\trateCT\tpopSize\n')

    i = 0
    sampleNr = 0

    with open(trees_path, 'r') as fp:
        for line in fp:

            line = line.rstrip('\n').rstrip('\r')
            tree = dendropy.Tree.get(string=line, schema="newick", tree_offset=0, taxon_namespace=taxa)

            sample = 0

            tree_file = open('{}{}.tree'.format(seq_path, i), "w")

            csv_file = open('{}{}.csv'.format(seq_path, i), "w")
            csv_file.write("sample\tposterior\tlikelihood\tprior\tCoalescentConstant\trootHeight\tfreqParameter.1\tfreqParameter.2\tfreqParameter.3\tfreqParameter.4")
            csv_file.write("\trateAC\trateAG\trateAT\trateCG\trateCT\tpopSize\n")

            with open('{}{}.log'.format(seq_path, i), 'r') as fp2:
                for line2 in fp2:
                    line2 = line2.rstrip('\n').rstrip('\r')
                    if line2.startswith('lp__'):
                        header = line2.split(',')
                        header.pop(0)
                    elif not line2.startswith('#') and not line2.startswith('lp__'):
                        elements = line2.split(',')
                        elements.pop(0)
                        if len(elements) == len(header):

                            params = {}
                            for j in range(len(elements)):
                                params[header[j]] = float(elements[j])

                            upadate_tree(tree, params)

                            freqs = [float(params['freqs.{}'.format(pi + 1)]) for pi in range(4)]
                            rates = [float(params['rates.{}'.format(r + 1)]) for r in range(6)]

                            model = phylo.GTR(rates, freqs)

                            try:
                                likelihood = phylo.compute_likelihood(tree, dna, model)

                                prior_coalescent = phylo.constant_coalescent(
                                    [float(params['taus.{}'.format(idx + 1)]) for idx in range(len(dna) - 1)],
                                    float(params['theta']))

                                prior_rates = scipy.stats.dirichlet(rates, [1] * 6)
                                prior_freqs = scipy.stats.dirichlet(freqs, [1] * 4)
                            except ValueError, e:
                                likelihood = 0
                                prior_freqs = 0
                                prior_rates = 0
                                prior_coalescent = 0

                            prior = prior_coalescent + prior_rates + prior_freqs
                            posterior = likelihood + prior

                            rates = map(lambda x: x / rates[len(rates) - 1])
                            rates.pop()

                            l = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(sample, posterior, likelihood, prior, prior_coalescent, 0)
                            l += '\t'.join([str(x) for x in freqs])
                            l += '\t' + '\t'.join([str(x) for x in rates])
                            l += '\t{}\n'.format(params['theta'])
                            csv_file.write(l)

                            output_probs.write('{}\t{}\t'.format(sampleNr, i))
                            output_probs.write(l)

                            sample += 1
                            sampleNr += 1
                            tree_file.write(tree.as_string(schema="newick"))
            tree_file.close()
            sys.stderr.write('tree {}\n'.format(i))
            i += 1
    output_probs.close()


def create_trees(stem_file, tree_file):
    taxa = dendropy.TaxonNamespace()

    tree_template = stem_file + '{}.tree'
    log_template = stem_file + '{}.log'
    i = 0
    with open(tree_file, 'r') as fp:
        for line in fp:
            line = line.rstrip('\n').rstrip('\r')
            tree = dendropy.Tree.get(string=line, schema="newick", tree_offset=0, taxon_namespace=taxa)

            output = open(tree_template, "w")

            with open(log_template.format(i), 'r') as fp2:
                for line2 in fp2:
                    line2 = line2.rstrip('\n').rstrip('\r')
                    if line2.startswith('lp__'):
                        header = line2.split(',')
                        header.pop(0)
                    elif not line2.startswith('#') and not line2.startswith('lp__'):
                        elements = line2.split(',')
                        elements.pop(0)
                        if len(elements) == len(header):

                            params = {}
                            for j in range(len(elements)):
                                params[header[j]] = float(elements[j])

                            upadate_tree(tree, params)

                            output.write(tree.as_string(schema="newick"))
            output.close()
            sys.stderr.write('tree {}\n'.format(i))
            i += 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('alignment')
    parser.add_argument('trees')
    args = parser.parse_args()

    create_trees(args.alignment, args.trees)
    recompute_all(args.alignment, args.trees)
