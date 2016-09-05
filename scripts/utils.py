import dendropy
import sys
import re
from sets import Set


def convert_trees(input, output=None):
    handle = open(output, 'w') if output else sys.stdout
    with open(input, 'r') as fp:
        for line in fp:
            line = line.rstrip('\n').rstrip('\r')
            line = line[1:len(line)-1]
            splits = re.findall(r'\[[^\]]+\]', line)

            n = max([int(x) for x in re.findall(r'\d+', line)])
            taxon_namespace = dendropy.TaxonNamespace(['t' + str(i+1) for i in range(n)])
            tree = dendropy.Tree(taxon_namespace=taxon_namespace)

            nodes = {}

            for idx in range(n):
                n = dendropy.Node(edge_length=0.0)
                n.taxon = taxon_namespace.get_taxon('t' + str(idx+1))
                nodes[str(idx+1)] = n

            current_height = 0.0
            for s in splits:
                current_height += 1.0
                idxs = re.findall(r'\d+', s)

                s = Set([])
                for i in idxs:
                    s.add(nodes[i])

                node = dendropy.Node(edge_length=current_height)
                for ss in s:
                    node.add_child(ss)
                for i in idxs:
                    nodes[i] = node
            tree.seed_node = nodes['1']

            for nd in tree.postorder_node_iter():
                if nd.parent_node is None:
                    nd.edge.length = 0.0
                else:
                    nd.edge.length = nd.parent_node.edge.length - nd.edge.length

            handle.write(tree.as_string(schema="newick"))

    if handle is not sys.stdout:
        handle.close()


if __name__ == '__main__':
    convert_trees('/Users/mathieu/Downloads/tTauCurvature-master/all_trees_on_6_taxa.txt')
