import dendropy
import argparse

from subprocess import check_call, STDOUT
from dendropy.simulate import treesim

parser = argparse.ArgumentParser()
parser.add_argument('tree')
parser.add_argument('alignment')
args = parser.parse_args()

taxa = dendropy.TaxonNamespace(["t1", "t2", "t3", "t4", "t5", "t6"])
tree = treesim.pure_kingman_tree(taxon_namespace=taxa, pop_size=0.2)
print(tree.as_string(schema="newick"))
print(tree.as_ascii_plot())

dest = open(args.tree, "w")
dest.write(tree.as_string(schema="newick"))
dest.close()

freqs = 0.33, 0.25, 0.16, 0.26
#rates = 0.04, 0.4, 0.02, 0.02, 0.5, 0.01
rates = 4, 40, 2, 2, 50

cmd = ['/Users/mathieu/CProjects/phyloc/phyloc/Release/simultron', '-i', args.tree, '-o', args.alignment, '-F', 'fasta', '-m', 'GTR', '-l', '2000']
cmd.extend(['-f', ','.join([str(x) for x in freqs])])
cmd.extend(['-r', ','.join([str(x) for x in rates])])
print(cmd)
check_call(cmd)