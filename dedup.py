#!/usr/bin/python
#-*- coding: UTF-8 -*-

import os
import sys
from optparse import OptionParser
import pysam
import time

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
    return parser.parse_args()


def str_compare(forward_list1, forward_list2, reverse_list1, reverse_list2, phredQ=20):
    # I believe this is the right comparing. but keep the same output with origin version first.
    # if C_EXT: return c_ext.str_compare(
    #     len(forward_list1[2]), forward_list1[2], forward_list1[3].tolist(),
    #     len(forward_list2[2]), forward_list2[2], forward_list2[3].tolist(),
    #     phredQ) and c_ext.str_compare(
    #     len(reverse_list1[1]), reverse_list1[1], reverse_list1[2].tolist(),
    #     len(reverse_list2[1]), reverse_list2[1], reverse_list2[2].tolist(),
    #     phredQ)

    seq_1 = forward_list1[2] + reverse_list1[1]
    phred_1 = forward_list1[3].tolist()
    tmp_2 = reverse_list1[2].tolist()
    phred_1.extend(tmp_2)
    seq_2 = forward_list2[2] + reverse_list2[1]
    phred_2 = forward_list2[3].tolist()
    tmp_2 = reverse_list2[2].tolist()
    phred_2.extend(tmp_2)

    if C_EXT: return c_ext.str_compare(
        len(seq_1), seq_1, phred_1,
        len(seq_2), seq_2, phred_2,
        phredQ)

    if len(seq_1) == len(seq_2):
        if seq_1 == seq_2:
            return True
        else:
            index = [0]
            count = 0
            for i in range(len(seq_1)):
                if phred_1[i] < phredQ or phred_2[i] < phredQ:
                    if i == 0:
                        if seq_1[i] != seq_2[i]:
                            return False
                    else:
                        index.append(i)
                        count += 1
                        if seq_1[(index[count - 1] + 1):index[count]] != seq_2[(index[count - 1] + 1):index[count]]:
                            return False
                else:
                    if seq_1[i] == seq_2[i]:
                        continue
                    else:
                        return False
            return True
    else:
        return False


def block_compare(forward_block_dict_list, reverse_block_dict, phredQ=20):
    dup = set()
    keep = [] # list of id of reads to keep
    for i in range(len(forward_block_dict_list)):
        if forward_block_dict_list[i][1] not in dup:
            keep.append(forward_block_dict_list[i][0])
            if forward_block_dict_list[i][1] not in reverse_block_dict:
                continue
        for j in range((i + 1), len(forward_block_dict_list)):
            if forward_block_dict_list[j][1] not in dup:
                if forward_block_dict_list[j][1] in reverse_block_dict:
                    if str_compare(forward_block_dict_list[i], forward_block_dict_list[j], reverse_block_dict[forward_block_dict_list[i][1]], reverse_block_dict[forward_block_dict_list[j][1]], phredQ):
                        dup.add(forward_block_dict_list[j][1])
    return keep, dup


def dedup(in_file, phred, out_file):
    f = pysam.AlignmentFile(in_file, "rb")
    forward_block_dict = {}
    reverse_block_dict = {}
    chrom_flag = False
    chrom_tmp = "chrM"
    dup = set() # duplicated pair names
    keep = [] # bitmap, whether the Nth read should be kept
    for rid, frag in enumerate(f):
        keep.append(False)
        if frag.is_proper_pair:
            if frag.template_length > 0:
                if frag.reference_id not in forward_block_dict:
                    forward_block_dict[frag.reference_id] = {}
                    forward_block_dict[frag.reference_id][frag.reference_start] = [
                        [rid, frag.query_name, frag.query_sequence, frag.query_qualities]]
                    if not chrom_flag:
                        chrom_tmp = frag.reference_id
                        chrom_flag = True
                        continue
                    for i in forward_block_dict[chrom_tmp].keys():
                        if len(forward_block_dict[chrom_tmp][i]) >= 2:
                            block_keep, block_dup = block_compare(forward_block_dict[chrom_tmp][i],
                                                      reverse_block_dict[chrom_tmp], phred)
                            dup.update(block_dup)
                            for i in block_keep:
                                keep[i] = True
                        else:
                            keep[forward_block_dict[chrom_tmp][i][0][0]] = True
                    forward_block_dict[chrom_tmp].clear()
                    for i in reverse_block_dict[chrom_tmp].keys():
                        if i not in dup:
                            keep[reverse_block_dict[chrom_tmp][i][0]] = True
                    reverse_block_dict[chrom_tmp].clear()
                    chrom_tmp = frag.reference_id
                elif frag.reference_start not in forward_block_dict[frag.reference_id]:
                    forward_block_dict[frag.reference_id][
                        frag.reference_start] = {}
                    forward_block_dict[frag.reference_id][frag.reference_start] = [
                        (rid, frag.query_name, frag.query_sequence, frag.query_qualities)]
                else:
                    forward_block_dict[frag.reference_id][frag.reference_start].append(
                        (rid, frag.query_name, frag.query_sequence, frag.query_qualities))
            elif frag.template_length < 0:
                if frag.reference_id not in reverse_block_dict:
                    reverse_block_dict[frag.reference_id] = {}
                reverse_block_dict[frag.reference_id][frag.query_name] = rid, frag.query_sequence, frag.query_qualities
            else:
                keep[rid] = True
        else:
            keep[rid] = True
# 最后一次reference_id变换后的数据处理 开始
    if chrom_tmp in forward_block_dict:  # chrom_tmp可能不存在
        for i in forward_block_dict[chrom_tmp].keys():
            if len(forward_block_dict[chrom_tmp][i]) >= 2:
                block_keep, block_dup = block_compare(forward_block_dict[chrom_tmp][
                                          i], reverse_block_dict[chrom_tmp], phred)
                dup.update(block_dup)
                for i in block_keep:
                    keep[i] = True
            else:
                keep[forward_block_dict[chrom_tmp][i][0][0]] = True
        forward_block_dict[chrom_tmp].clear()
        for i in reverse_block_dict[chrom_tmp].keys():
            if i not in dup:
                keep[reverse_block_dict[chrom_tmp][i][0]] = True
        reverse_block_dict[chrom_tmp].clear()
# 最后一次reference_id变换后的数据处理 结束
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
    dedup(o.input1, o.phred, o.output)
    time2 = time.time()
    print 'Time used: ' + str(time2 - time1)
