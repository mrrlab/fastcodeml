#!/usr/bin/env python
import sys
import argparse
import subprocess
import csv
import os


def convert_params(record):
    res = []
    for par in ('w0', 'k', 'p0', 'p1', 'w2'):
        if record[par] == 'NA':
            assert par == 'w2'
            continue
        res.extend(('--init-param',
                    '%s=%s' % (par, record[par])))
    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="continue fastcodeml from a specific point")
    parser.add_argument('parfile', help='tab-separated file with starting points')
    parser.add_argument('treefile', help='tree structure file name')
    parser.add_argument('seqfile', help='sequence data filename')
    parser.add_argument('outfile', help='main result file name')
    parser.add_argument('extra', nargs='*',
                        help='extra arguments to pass to the executable (use -- '
                        'before the first argument)')

    parser.add_argument('--exe', '-x', default='fast',
                        help='fastcodeml executable')


    args = parser.parse_args()

    reader = csv.DictReader(open(args.parfile), delimiter='\t')
    for record in reader:
        fargs = [args.exe, '--branch', record['branch'],
                 '--only-hyp', record['hyp']]
        init = convert_params(record)
        fargs.extend(init)
        fargs.extend(('--output',
                      '%s.%s.%s.fcml' % (args.outfile,
                                         record['branch'],
                                         record['hyp'])))

        fargs.extend(args.extra)
        fargs.extend((args.treefile, args.seqfile))
        subprocess.call(fargs)

