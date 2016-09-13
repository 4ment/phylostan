#!/usr/bin/env python
import utils
import os
import subprocess
import argparse

from subprocess import check_call

parser = argparse.ArgumentParser()
parser.add_argument('alignment')
parser.add_argument('trees')
args = parser.parse_args()

dir, file = os.path.split(args.alignment)

trees_path = os.path.join(dir, 'output.trees')
temp_tree_path = os.path.join(dir, 'temp.tree')

utils.convert_trees(args.trees, trees_path)

my_path = os.path.split(os.path.realpath(__file__))[0]

i = 0

with open(trees_path, 'r') as f:
    for line in f:
        line = line.rstrip('\n').rstrip('\r')
        with open(temp_tree_path, 'w') as f2:
            f2.write(line)

        output = args.alignment + str(i) + '.log'

        cmd = ['python', os.path.join(my_path, 'phystan.py'), '-t', temp_tree_path, '-i', args.alignment, '-o', output]
        print(' '.join(cmd))

        # Stan can fail so we just restart it until success
        done = False
        while not done:
            done = True
            try:
                check_call(cmd)
            except subprocess.CalledProcessError, e:
                done = False
        i += 1
