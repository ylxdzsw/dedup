#!/usr/bin/python
#-*- coding: UTF-8 -*-

import os
import sys
from optparse import OptionParser
import pysam
import time
from multiprocessing import Pool

# try load C extensions
try:
    from c_ext import lib as c_ext
    C_EXT = True
except ImportError:
    C_EXT = False

def parse_command():
    usage = "usage: ./dedup.py -1 /haplox/rawout/150626/150702/liliang_cfDNA-009_sort.bam -o /haplox/users/huang/mypy/liliang_cfDNA-009_sort_test.bam > liliang_cfDNA-009_sort_test_len.dedup"
    version = "%prog 1.0"
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-1", "--input1", dest="input1", help="the bam file")
    parser.add_option("-q", "--phredQ", dest="phred",
                      help="Phred Q", type="int", default="20")
    parser.add_option("-o", "--output", dest="output",
                      default="/haplox/users/fastq", help="output file")
    parser.add_option("-p", "--process", dest="process",
                      help="number of worker processes to start", type="int", default="0")
    return parser.parse_args()


def str_compare(forward1, forward2, reverse1, reverse2, phredQ=20):
    # I believe this is the right comparing. but keep the same output with origin version first.
    # if C_EXT: return c_ext.str_compare(
    #     len(forward1[2]), forward1[2], forward1[3].tolist(),
    #     len(forward2[2]), forward2[2], forward2[3].tolist(),
    #     phredQ) and c_ext.str_compare(
    #     len(reverse1[1]), reverse1[1], reverse1[2].tolist(),
    #     len(reverse2[1]), reverse2[1], reverse2[2].tolist(),
    #     phredQ)

    seq1   = forward1[2] + reverse1[1]
    phred1 = forward1[3].tolist() + reverse1[2].tolist()
    seq2   = forward2[2] + reverse2[1]
    phred2 = forward2[3].tolist() + reverse2[2].tolist()

    if C_EXT: return c_ext.str_compare(
        len(seq1), seq1, phred1,
        len(seq2), seq2, phred2,
        phredQ)

    if len(seq1) == len(seq2):
        if seq1 == seq2:
            return True
        else:
            index = [0]
            count = 0
            for i in range(len(seq1)):
                if phred1[i] < phredQ or phred2[i] < phredQ:
                    if i == 0:
                        if seq1[i] != seq2[i]:
                            return False
                    else:
                        index.append(i)
                        count += 1
                        if seq1[(index[count - 1] + 1):index[count]] != seq2[(index[count - 1] + 1):index[count]]:
                            return False
                else:
                    if seq1[i] == seq2[i]:
                        continue
                    else:
                        return False
            return True
    else:
        return False


def block_compare(forward_block, reverse_block, phredQ=20):
    dup = set()
    keep = [] # list of id of reads to keep
    for i in range(len(forward_block)):
        if forward_block[i][1] not in dup:
            keep.append(forward_block[i][0])
            if forward_block[i][1] not in reverse_block:
                continue
        for j in range((i + 1), len(forward_block)):
            if forward_block[j][1] not in dup:
                if forward_block[j][1] in reverse_block:
                    if str_compare(forward_block[i], forward_block[j], reverse_block[forward_block[i][1]], reverse_block[forward_block[j][1]], phredQ):
                        dup.add(forward_block[j][1])
    return keep, dup


def blocks(f, keep):
    forward_block = {} # {chrom: {startpos: [(rid, name, seq, qual)]}}
    reverse_block = {} # {chrom: {name: (rid, seq, qual)}}
    current_chrom = -1
    for rid, frag in enumerate(f):
        keep.append(False)

        if not frag.is_proper_pair:
            keep[rid] = True
            continue

        if frag.template_length > 0:
            if frag.reference_id not in forward_block:
                if current_chrom != -1:
                    yield forward_block.pop(current_chrom), reverse_block.pop(current_chrom)
                forward_block[frag.reference_id] = {}
                current_chrom = frag.reference_id

            if frag.reference_start not in forward_block[frag.reference_id]:
                forward_block[frag.reference_id][frag.reference_start] = [
                    (rid, frag.query_name, frag.query_sequence, frag.query_qualities)]
            else:
                forward_block[frag.reference_id][frag.reference_start].append(
                    (rid, frag.query_name, frag.query_sequence, frag.query_qualities))
        elif frag.template_length < 0:
            if frag.reference_id not in reverse_block:
                reverse_block[frag.reference_id] = {}
            reverse_block[frag.reference_id][frag.query_name] = rid, frag.query_sequence, frag.query_qualities
        else:
            keep[rid] = True

    # yield last one
    if current_chrom in forward_block:
        yield forward_block[current_chrom], reverse_block[current_chrom]


def dedup_block(forward_block, reverse_block, phred):
    dup, keep = set(), []

    for i in forward_block.keys():
        if len(forward_block[i]) >= 2:
            block_keep, block_dup = block_compare(forward_block[i], reverse_block, phred)
            dup.update(block_dup)
            keep += block_keep
        else:
            keep.append(forward_block[i][0][0])

    for i in reverse_block.keys():
        if i not in dup:
            keep.append(reverse_block[i][0])

    return keep


def dedup(in_file, phred, out_file, pool):
    f = pysam.AlignmentFile(in_file, "rb")
    keep = [] # bitmap, whether the Nth read should be kept
    task = [] # hold references to async tasks

    for forward_block, reverse_block in blocks(f, keep):
        if pool:
            task.append(pool.apply_async(dedup_block, (forward_block, reverse_block, phred)))
        else:
            for i in dedup_block(forward_block, reverse_block, phred):
                keep[i] = True

    if pool:
        for t in task:
            for i in t.get():
                keep[i] = True

    f.close()
    f = pysam.AlignmentFile(in_file, "rb")
    f_w = pysam.AlignmentFile(out_file, "wb", template=f)
    for rid, frag in enumerate(f):
        if keep[rid]:
            f_w.write(frag)
    f_w.close()
    print "total=", rid,  "---after dedup---", sum(keep)


if __name__ == "__main__":
    o, args = parse_command()
    if o.input1 == None or o.output == None:
        print "see -h for help"
        sys.exit(-1)
    time1 = time.time()
    pool = Pool(processes=o.process) if o.process > 0 else None
    dedup(o.input1, o.phred, o.output, pool)
    time2 = time.time()
    print 'Time used: ' + str(time2 - time1)
