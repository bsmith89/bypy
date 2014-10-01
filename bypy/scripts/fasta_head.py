#!/usr/bin/env python
# Cut and return the top sequences from a sequence file
#
# TODO: Fix the Broken Pipe Error that I get every time I try to Pipe
# this command with something else.

import optparse
import sys
import Bio.SeqIO

TAIL_TO_OUT_FORMAT_MAPPING = {'fasta-tail':'fasta', 'fasta':'fasta'}
DEFAULT_N = 10
DEFAULT_FORMAT = 'fasta'

def write_head(seq_file, num, frmt):
    recs_iter = Bio.SeqIO.parse(seq_file, frmt)
    recs_out = []
    for i, rec in enumerate(recs_iter):
        if i == num:
            break
        else:
            recs_out += [rec]
    Bio.SeqIO.write(recs_out, sys.stdout, frmt)
        

def write_tail(seq_file, num, frmt):
    out_frmt = TAIL_TO_OUT_FORMAT_MAPPING[frmt]
    recs_iter = Bio.SeqIO.parse(seq_file, frmt)
    recs_out = []
    for i, rec in enumerate(recs_iter):
        if i == num:
            break
        else:
            recs_out += [rec]
    Bio.SeqIO.write(recs_out, sys.stdout, out_frmt)
    


def main():
    usage = "usage: %prog [options] file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-n", "--number", "--num", dest="num_seqs",
                      type='int', default=DEFAULT_N,
                      help="number of sequences to retrieve. \
                            DEFAULT: %d" % DEFAULT_N)
    parser.add_option("-t", "--tail", "--bottom-up", action="store_true",
                      dest="tail",
                      help="take from the bottom of the file instead of the top. \
                            DEFAULT: Take from the top")
    parser.add_option("-f", "--format", "--file-format", "--file-type",
                      dest="file_format", default=DEFAULT_FORMAT,
                      help="the format of the sequence file. \
                            Can take any sequence string accepted by biopython. \
                            DEFAULT: %s" % DEFAULT_FORMAT)
    (opts, args) = parser.parse_args()
    opts.file_format
    if not opts.tail:
        if len(args) > 0:
            for seq_file_path in args:
                write_head(seq_file_path, opts.num_seqs, opts.file_format)
        else:
            write_head(sys.stdin, opts.num_seqs, opts.file_format)
    else:
        if len(args) > 0:
            for seq_file_path in args:
                write_tail(seq_file_path, opts.num_seqs, opts.file_format)
        else:
            write_tail(sys.stdin, opts.num_seqs, opts.file_format)


if __name__ == "__main__":
    main()
