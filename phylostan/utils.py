import csv

import numpy as np


def ratios_root_height_from_branch_lengths(tree, eps=1.0e-6):
    bounds = get_lowers(tree)
    heights = heights_from_branch_lengths(tree, eps)
    tip_count = len(tree.taxon_namespace)
    ratios = [None] * (tip_count - 2)
    for node in tree.preorder_node_iter(
        lambda n: n != tree.seed_node and n.is_internal()
    ):
        ratios[node.index - 1 - tip_count] = (
            heights[node.index - 1] - bounds[node.index - 1]
        ) / (heights[node.parent_node.index - 1] - bounds[node.index - 1])
    return np.array(ratios), np.array(heights[-1:])


def heights_from_branch_lengths(tree, eps=1.0e-6):
    heights = [None] * (2 * len(tree.taxon_namespace) - 1)
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            heights[node.index - 1] = node.date
        else:
            heights[node.index - 1] = max(
                [
                    heights[c.index - 1] + max(eps, c.edge_length)
                    for c in node.child_node_iter()
                ]
            )
    return heights


def setup_dates(tree, dates=None, heterochronous=False):
    # if a date file is provided then it is heterochronous
    if dates:
        heterochronous = True

    # parse dates
    if heterochronous:
        dates = {}
        if dates:
            with dates as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    dates[row['name']] = float(row['date'].strip())
        else:
            for node in tree.leaf_node_iter():
                dates[str(node.taxon)] = float(str(node.taxon).split('_')[-1][:-1])

        max_date = max(dates.values())
        min_date = min(dates.values())

        # time starts at 0
        if min_date == 0:
            for node in tree.leaf_node_iter():
                node.date = dates[str(node.taxon)]
            oldest = max_date
        # time is a year
        else:
            for node in tree.leaf_node_iter():
                node.date = max_date - dates[str(node.taxon)]
            oldest = max_date - min_date
    else:
        for node in tree.postorder_node_iter():
            node.date = 0.0
        oldest = None

    return oldest


def setup_indexes(tree):
    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = len(tree.taxon_namespace) + 1
    taxa_dict = {taxon.label: idx for idx, taxon in enumerate(tree.taxon_namespace)}

    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            s += 1
        else:
            node.index = taxa_dict[node.taxon.label] + 1


def get_peeling_order(tree):
    peeling = []
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            peeling.append([x.index for x in node.child_node_iter()] + [node.index])
    return peeling


def get_preorder(tree):
    peeling = []
    peeling.append([tree.seed_node.index, 0])
    for node in tree.preorder_node_iter():
        if node.parent_node is not None:
            peeling.append([node.index, node.parent_node.index])
    return peeling


def get_lowers(tree):
    lowers = [0 for x in tree.postorder_node_iter()]
    ll = {}

    for node in tree.postorder_node_iter():
        if node.is_leaf():
            ll[node] = node.date
        else:
            ll[node] = max([ll[x] for x in node.child_node_iter()])

    for node in tree.preorder_node_iter():
        lowers[node.index - 1] = ll[node]
    return lowers


def get_dna_leaves_partials(alignment):
    tipdata = np.zeros((len(alignment), alignment.sequence_size, 4), dtype=np.int)
    dna_map = {
        'a': [1, 0, 0, 0],
        'c': [0, 1, 0, 0],
        'g': [0, 0, 1, 0],
        't': [0, 0, 0, 1],
    }

    s = 0
    for name in alignment:
        for i, c in enumerate(alignment[name].symbols_as_string()):
            tipdata[s][i][:] = dna_map.get(c.lower(), [1, 1, 1, 1])
        s += 1
    return tipdata


def get_dna_leaves_partials_compressed(alignment):
    weights = []
    keep = [True] * alignment.sequence_size

    patterns = {}
    indexes = {}
    for i in range(alignment.sequence_size):
        pat = ''
        for name in alignment:
            pat += str(alignment[name][i])

        if pat in patterns:
            keep[i] = False
            patterns[pat] += 1
        else:
            patterns[pat] = 1
            indexes[i] = pat

    for i in range(alignment.sequence_size):
        if keep[i]:
            weights.append(patterns[indexes[i]])

    tipdata = np.zeros((len(alignment), len(patterns.keys()), 4), dtype=np.int)

    dna_map = {
        'a': [1, 0, 0, 0],
        'c': [0, 1, 0, 0],
        'g': [0, 0, 1, 0],
        't': [0, 0, 0, 1],
    }

    s = 0
    for name in alignment:
        k = 0
        for i, c in enumerate(alignment[name].symbols_as_string()):
            if keep[i]:
                tipdata[s][k][:] = dna_map.get(c.lower(), [1, 1, 1, 1])
                k += 1
        s += 1
    return tipdata, weights


def to_nexus(node, fp):
    if not node.is_leaf():
        fp.write('(')
        for i, n in enumerate(node.child_node_iter()):
            to_nexus(n, fp)
            if i == 0:
                fp.write(',')
        fp.write(')')
    else:
        fp.write(str(node.index))
    if hasattr(node, 'date'):
        fp.write('[&height={}'.format(node.date))
        if hasattr(node, 'rate'):
            fp.write(',rate={}'.format(node.rate))
        fp.write(']')
    if node.parent_node is not None:
        fp.write(':{}'.format(node.edge_length))
    else:
        fp.write(';')


def convert_samples_to_nexus(tree, input, output, rate=None):
    taxaCount = len(tree.taxon_namespace)
    outp = open(output, 'w')
    outp.write('#NEXUS\nBegin trees;\nTranslate\n')
    outp.write(
        ',\n'.join(
            [
                str(i + 1) + ' ' + x.label.replace("'", '')
                for i, x in enumerate(tree.taxon_namespace)
            ]
        )
    )
    outp.write('\n;\n')

    time = False
    with open(input) as fp:
        for line in fp:
            if line.startswith('lp'):
                try:
                    line.split(',').index('blens.1')
                except ValueError:
                    time = True
                break

    count = 1
    if time:
        with open(input) as fp:
            for line in fp:
                if line.startswith('lp'):
                    header = line.split(',')
                    hindex = header.index('heights.1')
                    if rate is not None:
                        strict = False
                        try:
                            rindex = header.index('substrates.1')
                        except ValueError:
                            rindex = header.index('rate')
                            strict = True
                    else:
                        strict = True
                elif not line.startswith('#'):
                    l = line.split(',')
                    for n in tree.postorder_node_iter():
                        if not n.is_leaf():
                            n.date = float(l[hindex + n.index - taxaCount - 1])
                    if rate is not None:
                        if strict:
                            for n in tree.postorder_node_iter():
                                if n.parent_node is not None:
                                    n.rate = float(l[rindex]) if rate is None else rate
                        else:
                            for n in tree.postorder_node_iter():
                                if n.parent_node is not None:
                                    n.rate = float(l[rindex + n.index - 1])

                    for n in tree.postorder_node_iter():
                        if n.parent_node is not None:
                            n.edge_length = n.parent_node.date - n.date
                    outp.write('tree {} = '.format(count))
                    count += 1
                    to_nexus(tree.seed_node, outp)
                    outp.write('\n')
    else:
        with open(input) as fp:
            for line in fp:
                if line.startswith('lp'):
                    header = line.split(',')
                    bindex = header.index('blens.1')
                elif not line.startswith('#'):
                    l = line.split(',')
                    for n in tree.postorder_node_iter():
                        if n.parent_node is not None:
                            if bindex + n.index - 1 < len(l):
                                n.edge_length = float(l[bindex + n.index - 1])
                            else:
                                n.edge_length = 0.0
                    outp.write('tree {} = '.format(count))
                    count += 1
                    to_nexus(tree.seed_node, outp)
                    outp.write('\n')
    outp.write('END;')
    outp.close()


def descriptive_stats(d, alpha):
    median, low, high = np.quantile(d, (0.5, alpha / 2.0, 1.0 - alpha / 2.0))
    return np.mean(d), median, low, high


def parse_log(inputfile, alpha=0.05, tree=None):
    data = []
    GTR = ('AC', 'AG', 'AT', 'CG', 'CT', 'GC')
    frequencies = 'A', 'C', 'G', 'T'
    variables = {
        'wshape': 'Weibull (shape)',
        'pinv': 'Proportion invariant',
        'kappa': 'HKY (kappa)',
        'rate': 'Strict clock (rate)',
        'theta': 'Constant population size (theta)',
        'tau': 'GMRF precision (tau)',
        'netDiversificationRate': 'net diversification rate',
        'relativeExtinctionRate': 'relative extinction rate',
        'ucln_mean': 'UCLN mean',
        'ucln_stdev': 'UCLN stdev',
    }
    with open(inputfile) as fp:
        for line in fp:
            line = line.strip()
            if line.startswith('lp'):
                header = line.split(',')
                [data.append([]) for _ in range(len(header))]
            elif not line.startswith('#') and len(line) != 0:
                for idx, h in enumerate(line.split(',')):
                    data[idx].append(float(h))
        for var in header:
            # GTR
            if var == 'rates.1':
                print('GTR')
                for i in range(6):
                    d = data[header.index('rates.' + str(i + 1))]
                    mean, median, low, high = descriptive_stats(d, alpha)
                    print(
                        '  {} mean: {:.3E} {}% CI: ({:.3E},{:.3E})'.format(
                            GTR[i], mean, (1 - alpha) * 100, low, high
                        )
                    )
                for i in range(4):
                    d = data[header.index('freqs.' + str(i + 1))]
                    mean, median, low, high = descriptive_stats(d, alpha)
                    print(
                        '  {} mean: {:.4f} {}% CI: ({:.4f},{:.4f})'.format(
                            frequencies[i], mean, (1 - alpha) * 100, low, high
                        )
                    )
            elif var in variables:
                if 'substrates.1' not in header or var != 'rate':
                    d = data[header.index(var)]
                    mean, median, low, high = descriptive_stats(d, alpha)
                    print(
                        '{} mean: {} {}% CI: ({},{})'.format(
                            variables[var], mean, (1 - alpha) * 100, low, high
                        )
                    )

            # time tree
        if 'heights.1' in header:
            index_root_height = 0
            max_height = 0
            for idx, h in enumerate(header):
                if h.startswith('heights.') and data[idx][0] > max_height:
                    index_root_height = idx
                    max_height = data[idx][0]
            d = data[index_root_height]
            mean, median, low, high = descriptive_stats(d, alpha)
            print(
                'Root height mean: {} {}% CI: ({},{})'.format(
                    mean, (1 - alpha) * 100, low, high
                )
            )
            if tree is not None and 'substrates.1' in header:
                rindex = header.index('substrates.1')
                hindex = header.index('heights.1')
                taxaCount = len(tree.taxon_namespace)
                d = []
                rates = []
                variances = []
                for i in range(len(data[0])):
                    for n in tree.postorder_node_iter():
                        if not n.is_leaf():
                            n.date = data[hindex + n.index - taxaCount - 1][i]
                        if n.parent_node is not None:
                            n.rate = data[rindex + n.index - 1][i]
                            rates.append(n.rate)
                    sum_distance = 0
                    sum_time = 0
                    for n in tree.postorder_node_iter():
                        if n.parent_node is not None:
                            edge_length = n.parent_node.date - n.date
                            sum_distance += n.rate * edge_length
                            sum_time += edge_length
                    d.append(sum_distance / sum_time)
                    variances.append(np.var(rates))
                mean, median, low, high = descriptive_stats(d, alpha)
                print(
                    'Mean rate mean: {} {}% CI: ({},{})'.format(
                        mean, (1 - alpha) * 100, low, high
                    )
                )
                mean, median, low, high = descriptive_stats(variances, alpha)
                print(
                    'Variance rate mean: {} {}% CI: ({},{})'.format(
                        mean, (1 - alpha) * 100, low, high
                    )
                )
        else:
            indexes = []
            for idx, h in enumerate(header):
                if h.startswith('blens'):
                    indexes.append(idx)
            sums = []
            for row in range(len(data[0])):
                sum_blens = 0
                for idx in indexes:
                    sum_blens += data[idx][row]
                sums.append(sum_blens)
            mean, median, low, high = descriptive_stats(sums, alpha)
            print(
                'Tree length mean: {} {}% CI: ({},{})'.format(
                    mean, (1 - alpha) * 100, low, high
                )
            )
