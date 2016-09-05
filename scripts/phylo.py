import numpy
import math


class GTR:
    """
    GTR substitution model
    """
    def __init__(self, rates, pi):
        self.rates = rates
        self.pi = pi

        self.Q = numpy.zeros((4, 4))

        self.eigvals = numpy.zeros(4)
        self.eigvecs = numpy.zeros((4, 4))
        self.eigvecs_inv = numpy.zeros((4, 4))

    def get_pi(self, i):
        return self.pi[i]

    def p_t(self, t):
        return numpy.abs(numpy.dot(self.eigvecs, numpy.dot(numpy.diag(numpy.exp(self.eigvals * t)), self.eigvecs_inv)))

    def update(self):
        self.update_Q()
        self.update_eigen_system()

    def update_Q(self):
        self.Q[0, 1] = self.Q[1, 0] = self.rates[0]  # a
        self.Q[0, 2] = self.Q[2, 0] = self.rates[1]  # b
        self.Q[0, 3] = self.Q[3, 0] = self.rates[2]  # c

        self.Q[1, 2] = self.Q[2, 1] = self.rates[3]  # d
        self.Q[1, 3] = self.Q[3, 1] = self.rates[4]  # e

        # should sum to 1
        if len(self.rates) == 6:
            self.Q[2, 3] = self.Q[3, 2] = self.rates[5]  # f
        else:
            self.Q[2, 3] = self.Q[3, 2] = 1.0  # f

        for i in range(4):
            for j in xrange(i + 1, 4):
                temp = self.Q[i, j]
                self.Q[i, j] = temp * self.pi[j]
                self.Q[j, i] = temp * self.pi[i]

    def update_eigen_system(self):
        numpy.fill_diagonal(self.Q, 0.0)
        for i in range(4):
            self.Q[i, i] = -numpy.sum(self.Q[i])

        subst = 0.0
        for i in range(4):
            subst -= self.Q[i, i] * self.pi[i]
        self.Q /= subst

        eigvals, eigvecs = numpy.linalg.eig(self.Q)
        self.eigvals = eigvals
        self.eigvecs = eigvecs
        self.eigvecs_inv = numpy.linalg.inv(self.eigvecs)


def get_peeling_order(tree):
    peeling = []
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            peeling.append(
                [x.index for x in node.child_node_iter()] + [node.index])
    return peeling


def get_dna_leaves_partials(alignment):
    tipdata = numpy.zeros((len(alignment), alignment.sequence_size, 4), dtype=numpy.int)
    dna_map = {'a': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'g': [0, 0, 1, 0], 't': [0, 0, 0, 1]}

    s = 0
    for name in alignment:
        for i, c in enumerate(alignment[name].symbols_as_string()):
            tipdata[s][i][:] = dna_map.get(c.lower(), [1, 1, 1, 1])
        s += 1
    return tipdata


def initialize_dna_partials(alignment):
    partials = numpy.zeros((2 * len(alignment) - 1, alignment.sequence_size, 4))
    dna_map = {'a': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'g': [0, 0, 1, 0], 't': [0, 0, 0, 1]}

    s = 0
    for name in alignment:
        for i, c in enumerate(alignment[name].symbols_as_string()):
            partials[s][i][:] = dna_map.get(c.lower(), [1, 1, 1, 1])
        s += 1
    return partials


def setup_indexes(tree, dna):
    """
    Set up indexes for each node of tree.
    A leaf node index is its position in the alignment.
    :param tree:
    :param dna:
    :return:
    """
    sequence_count = len(dna)

    for node in tree.postorder_node_iter():
        node.index = -1
        node.annotations.add_bound_attribute("index")

    s = sequence_count
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            node.index = s
            s += 1
        else:
            for idx, name in enumerate(dna):
                if str(name) == str(node.taxon):
                    node.index = idx
                    break


def constant_coalescent(times, pop):
    """
    Calculate log probability of intercoalescent intervals and population size
    :param times: intercoalescent intervals
    :param pop: population size
    :return: log probability
    """
    p = 0
    l = len(times) + 1.0

    for time in times:
        temp = ((l * (l - 1.0)) / 2.0) / pop
        #p += math.log(temp) - (temp * time)
        p -= math.log(pop) + (temp * time)
        l -= 1.0
    return p


def traverse(node, matrices, partials, model):
    """
    Traverse a tree in post order and update partial likelihoods and transition rate matrices
    :param node: current node
    :param matrices: array of matrices
    :param partials: array of partials
    :param model: substitution model
    :return:
    """

    idx = node.index
    if node.parent_node is not None:
        matrices[idx] = model.p_t(node.edge.length)

    if node.is_internal():
        indexes = []
        for c in node.child_node_iter():
            indexes.append(c.index)
            traverse(c, matrices, partials, model)

        sequence_length = len(partials[0])
        for i in range(sequence_length):
            partials[idx, i] = numpy.dot(matrices[indexes[0]], partials[indexes[0], i]) * numpy.dot(
                matrices[indexes[1]], partials[indexes[1], i])


def compute_likelihood(tree, alignment, model):
    """
    Calculate the log probability of alignment given tree and substitution model
    :param tree: phylogenetic tree
    :param alignment: aligned sequences
    :param model: substitution model
    :return: log likelihood
    """

    setup_indexes(tree, alignment)

    model.update()

    partials = initialize_dna_partials(alignment)

    matrices = numpy.zeros((2 * len(alignment) - 2, 4, 4))

    traverse(tree.seed_node, matrices, partials, model)

    lnl = 0
    for i in range(alignment.sequence_size):
        p = 0.0
        for j in range(4):
            p += partials[2 * len(alignment) - 2][i][j] * model.get_pi(j)
        lnl += math.log(p)
        # try:
        #     lnl += math.log(p)
        # except ValueError, e:
        #     sys.stderr.write(e)
        #     pass
    return lnl


