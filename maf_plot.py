#!/usr/bin/env python
'''
  plot maf files
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

plt.rc('savefig', dpi=300)
FIGSIZE=(16,12)

def main(prefix, out):
  # mutation frequencies
  first = True
  sample = []
  silent = []
  nonsilent = []
  for row in csv.reader(open('{}.sample_mutation_rate.tsv'.format(prefix), 'r'), delimiter='\t'):
    if first:
      first = False
      continue
    sample.append(row[0])
    silent.append(float(row[1]))
    nonsilent.append(float(row[2]))

  # sort from high to low
  indexes = sorted(range(len(silent)), key=lambda k: -silent[k])

  sample_sorted = [sample[x] for x in indexes]
  silent_sorted = [silent[x] for x in indexes]
  nonsilent_sorted = [nonsilent[x] for x in indexes]

  fig = plt.figure()
  
  ax = fig.add_subplot(111)
  ax.scatter(sample_sorted, silent_sorted, label='Silent mutations', s=4, c='#801010')
  ax.scatter(sample_sorted, nonsilent_sorted, label='Non-silent mutations', s=4, c='#303030')
  ax.set_xlabel('Samples')
  ax.set_ylabel('Mutation rate (per Mb)')
  ax.set_title('Mutation frequency')
  ax.set_yscale('log')
  ax.set_xticklabels([])
  ax.legend()

  fig.savefig('{}.sample_mutation_rate.png'.format(out))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='plot maf files')
  parser.add_argument('--prefix', default='stats', help='input filename prefix')
  parser.add_argument('--out', default='stats', help='filename prefix')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.prefix, args.out)
