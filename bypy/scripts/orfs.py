#!/usr/bin/env python

'''
File: extract_orfs.py
Author: Byron J Smith
Description: A multithreaded ORF translator from a nucleotide fasta file.

'''

import os
import sys
import optparse
import time
from Queue import Empty #, Queue
from multiprocessing import Process, cpu_count, Lock, Queue
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

_DEFAULT_IN_FMT = 'fasta'
_DEFAULT_OUT_FMT = 'fasta'
FILE_EXTENSION_MAP = {'.fastq': 'fastq',
                      '.fq': 'fastq',
                      '.fasta': 'fasta',
                      '.fa': 'fasta',
                      '.fn': 'fasta',
                      '.nucl': 'fasta',
                      '.aa': 'fasta'}
_DEFAULT_MIN_ORF_LEN = 20
_DEFAULT_TRANSLATE = False
_CODON_TABLE = 11
_TIMEOUT = 1
_DEFAULT_PROC_THREADS = 1
_WRITE_BUFF_SIZE = 100
##_MAX_QUEUE_SIZE = 10 * _WRITE_BUFF_SIZE
_MAX_QUEUE_SIZE = 10000


def calc_orfs(rec, min_orf_len=_DEFAULT_MIN_ORF_LEN,
              translate=_DEFAULT_TRANSLATE):
    if len(rec) / 3 < min_orf_len:
        return
    for strand, nucl_seq in [(+1, rec.seq),
                             (-1, rec.seq.reverse_complement())]:
        for frame in range(3):
            in_frame_nucl_seq = nucl_seq[frame:]
            start_i = 0
            for prot in in_frame_nucl_seq.translate(_CODON_TABLE).split("*"):
                end_i = start_i + (len(prot) * 3) + 3
                if end_i > len(in_frame_nucl_seq) - 3:
                    # We might just not have a stop codon, so when we're on
                    # the last possible codon we should not assume that the
                    # last codon exists and that it is a stop codon, the way
                    # we do for every other ORF.
                    end_i -= 3
                if len(prot) >= min_orf_len:
                    orf_nucl_seq = in_frame_nucl_seq[start_i:end_i]
                    if strand == +1:
                        first_nucl = start_i + frame
                        last_nucl = end_i + frame
                    elif strand == -1:
                        first_nucl = len(rec.seq) - frame - start_i
                        last_nucl = len(rec.seq) - frame - end_i
                    id = rec.id + "(%d)%d:%d" % (strand * (frame + 1),
                                                 first_nucl, last_nucl)
                    if translate:
                        orf_seq = prot
                    else:
                        orf_seq = orf_nucl_seq
                    orf_rec = SeqRecord(orf_seq, id, name='', description='')
                    yield orf_rec
                start_i = end_i


def read_recs(queue, path, file_format, quit_lock):
    if path is None:
        handle = sys.stdin
    else:
        handle = open(path)
    rec_iter = SeqIO.parse(path, file_format)
    for rec in rec_iter:
        queue.put(rec, block=True, timeout=None)
    quit_lock.release()
    if path is not None:
        handle.close()


def process_orfs(in_queue, out_queue, quit_lock,
                 min_orf_len=_DEFAULT_MIN_ORF_LEN,
                 translate=_DEFAULT_TRANSLATE):
    while True:
        try:
            rec = in_queue.get(block=True, timeout=_TIMEOUT)
        except Empty:
            if quit_lock.acquire(block=False):
                quit_lock.release()
                break
            else:
                continue
        else:
            for orf_rec in calc_orfs(rec, min_orf_len=min_orf_len,
                                     translate=translate):
                out_queue.put(orf_rec)


def write_orfs(queue, path, file_format, write_lock, quit_lock,
               buffer_size=_WRITE_BUFF_SIZE):
    if path is None:
        handle = sys.stdout  # It's not *really* the path, but it works.
    else:
        handle = open(path, 'w')
    rec_list = []
    while True:
        if len(rec_list) >= buffer_size:
            write_lock.acquire()
            SeqIO.write(rec_list, handle, file_format)
            write_lock.release()
            rec_list = []
        try:
            rec = queue.get(block=True, timeout=_TIMEOUT)
        except Empty:
            if quit_lock.acquire(block=False):
                quit_lock.release()
                write_lock.acquire()
                SeqIO.write(rec_list, handle, file_format)
                write_lock.release()
                if path is not None:
                    handle.close()
                break
            else:
                continue
        else:
            rec_list.append(rec)


def _determine_format_by_extension(path):
    if path is None:
        return None
    try:
        frmt = FILE_EXTENSION_MAP[os.path.splitext(path)[1]]
    except KeyError:
        frmt = None
    else:
        return frmt


def _determine_paths_formats(opts, args):
    if len(args) == 0:
        in_path = None
        out_path = None
    elif len(args) == 1:
        in_path = args[0]
        out_path = None
    elif len(args) == 2:
        in_path = args[0]
        out_path = args[1]
    else:
        raise ValueError("%d arguments passed. Accepts maximum of 2" %
                         len(args))
    in_format = opts.in_format
    if in_format is None:
        in_format = _determine_format_by_extension(in_path)
    if in_format is None:
        in_format = _DEFAULT_IN_FMT
    out_format = opts.out_format
    if out_format is None:
        out_format = _determine_format_by_extension(out_path)
    if out_format is None:
        out_format = _DEFAULT_OUT_FMT
    return in_path, in_format, out_path, out_format


def main():
    usage = "usage: %prog [options] [INFILE] [OUTFILE]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-l", "--min-length", "--min-orf-len",
                      dest="min_orf_len", default=_DEFAULT_MIN_ORF_LEN,
                      type='int',
                      help="the minimum ORF length. \
                            DEFAULT: %d" % _DEFAULT_MIN_ORF_LEN)
    parser.add_option("-p", "--cpus", dest="num_proc_threads",
                      default=_DEFAULT_PROC_THREADS, type='int',
                      help="number of processing threads to use. \
                            DEFAULT: %d" % _DEFAULT_PROC_THREADS)
    parser.add_option("-t", "--translate", dest="translate",
                      default=_DEFAULT_TRANSLATE, action='store_true',
                      help="translate the ORF. \
                            DEFAULT: %s" % _DEFAULT_TRANSLATE)
    parser.add_option("-T", "--do-not-translate", dest="translate",
                      default=_DEFAULT_TRANSLATE, action='store_false',
                      help="don't translate the ORF. \
                            DEFAULT: %s" % (not _DEFAULT_TRANSLATE))
    parser.add_option("-f", "--in-format", dest="in_format",
                      default=None,
                      help="format of the INFILE.\
                            Should be a valid Biopython format string.\
                            If not detected from extension, DEFAULT: %s" %
                           _DEFAULT_IN_FMT)
    parser.add_option("--out-format", dest="out_format",
                      default=None,
                      help="format of the OUTFILE.\
                            Should be a valid Biopython format string.\
                            If not detected from extension, DEFAULT: %s" %
                           _DEFAULT_OUT_FMT)
    parser.add_option("-n", "--num-out-files", dest="num_out_files",
                      default=1, type='int',
                      help="The number of out files to produce.\
                            DEFAULT: 1")
    opts, args = parser.parse_args()

    in_path, in_format, out_path, out_format = \
        _determine_paths_formats(opts, args)
    if opts.num_out_files > 1:
        base, extension = os.path.splitext(out_path)
        out_paths = [''.join([base, ".%02d" % i, extension])
                     for i in range(opts.num_out_files)]
    else:
        out_paths = [out_path]

    read_queue = Queue(maxsize=_MAX_QUEUE_SIZE)
    write_queue = Queue(maxsize=_MAX_QUEUE_SIZE)
    quit_lock = Lock()
    write_locks = [Lock() for i in range(opts.num_out_files)]
    read_recs_proc = Process(target=read_recs, args=(read_queue, in_path,
                                                    in_format, quit_lock))
    process_orfs_procs = [Process(target=process_orfs,
                                  args=(read_queue, write_queue,
                                        quit_lock, opts.min_orf_len,
                                        opts.translate))
                          for i in range(opts.num_proc_threads)]
    write_orfs_procs = [Process(target=write_orfs, args=(write_queue,
                                                        out_paths[i],
                                                        out_format,
                                                        write_locks[i],
                                                        quit_lock))
                        for i in range(opts.num_out_files)]
    quit_lock.acquire()
    read_recs_proc.start()
    for proc in process_orfs_procs:
        proc.start()
    for proc in write_orfs_procs:
        proc.start()
    all_procs = [read_recs_proc] + process_orfs_procs + write_orfs_procs

#    report_time = 0
#    while True:
#        print "Read Queue:",
#        print read_queue.qsize()
#        print "Write Queue:",
#        print write_queue.qsize()
#        time.sleep(2)
#        report_time += 2
#        if report_time == 60:
#            break

    for proc in all_procs:
        proc.join()

if __name__ == "__main__":
    main()
