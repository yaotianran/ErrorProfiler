#!/usr/bin/env python3
# v1.0a
# 将queryname sam文件按照模板长度分组
# 输出 filtered read count / total read count, mismatch base count / total base count
# filter_overlap_length_group_2.py <queryname sorted>.sam

import os
import os.path as path
import sys
import re
import pysam

SOFTCLIP_COMPILE = re.compile(r'\d+(?=S)')
SOFTCLIP_LEFT_COMPILE = re.compile(r'^\d+(?=S)')
SOFTCLIP_RIGHT_COMPILE = re.compile(r'\d+(?=S$)')
MATCH_COMPILE = re.compile(r'\d+(?=M)')
MD_COMPILE = re.compile(r'(?<=MD:Z:)\w+')
MD_SPLIT_COMPILE = re.compile(r'\d+|[A-Z]')

# 遍历一个queryname sorted sam文件，返回一对正常模式（Unambiguous scenario）的read string
# for forward_read, reverse_read in iter_queryname_sam(sam_file):
#     do samething about forward_read and reverse_read
#     line_str, 'header'

# 一对read string有两种配对模式，取决于两个read的map区域的绝对坐标的排序(S表示softclip)：
# 1. 正常模式（Unambiguous scenario）
# forward:   SSSS>>========================>>SSSSS
# reverse:                            SSSS<<====================<<SSSS
# 在这种模式下，两个read的TLEN的绝对值一样，forward read为正，reverse read为负

# 2. 交错模式(Ambiguous scenario)
# forward:                               SSS>>======================>>SSSS
# reverse:   SSSS<<========================<<SSSS
# 在这种模式下，两个read的TLEN的绝对值不一样
def iter_queryname_sam(sam_file: str, min_MAPQ = 0) -> tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
   '''
   提取两个read的overlap坐标和序列等信息

   Parameters:
      **sam_file**: str
         read pair的read string

   Returns:
      一对pysam.AlignedSegment forward and reverse
   '''
   samfile_str = path.realpath(path.expanduser(sam_file))
   bam_af = pysam.AlignmentFile(samfile_str)
   while True:
      try:

         segment1 = next(bam_af)
         segment2 = next(bam_af)

         while segment1.query_name != segment2.query_name:
            segment1 = segment2
            segment2 = next(bam_af)

         if segment1.mapping_quality < min_MAPQ or segment2.mapping_quality < min_MAPQ:
            continue

         if abs(segment1.template_length) != abs(segment1.template_length):
            continue


         if 'D' in segment1.cigarstring or 'I' in segment1.cigarstring or 'D' in segment2.cigarstring or 'I' in segment2.cigarstring: # cigar内有插入缺失
            continue

         if segment1.flag == 99 and segment2.flag == 147 or segment1.flag == 163 and segment2.flag == 83:  # read1为forward，read2为reverse
            yield segment1, segment2
            continue

         if segment1.flag == 147 and segment2.flag == 99 or segment1.flag == 83 and segment2.flag == 163: # read1为reverse，read2为forward
            yield segment2, segment1
            continue

      except StopIteration as ex:
         print('iter_queryname_sam: 已读取全部文件')
         return

      except Exception as ex:
         print('iter_queryname_sam:', segment1.query_name, '\n', segment2.query_name, 'not properly paired')
         print('iter_queryname_sam:', ex)
         print('skip')
         continue

   return




def get_base_PCR_cycle(segment_forward: pysam.AlignedSegment, segment_reverse: pysam.AlignedSegment) -> tuple[list, list]:
   '''
   得到一对reads中所有碱基在模板中的位置（PCR cycle）

   Parameters:
      **segment_forward， segment_reverse**: 一对pysam.AlignedSegment
         forward and reverse segments

   Returns: tuple[list, list]
      每个list由数值构成，第一个list为 matched bases的cycle (0-based)，第二个为mismatch bases的cycle (0-based)
   '''

   matched_list = []
   mismatch_list = []
   if segment_forward.flag == 99 and segment_reverse.flag == 147:  #  F1R2

      first_pair_list = segment_forward.get_aligned_pairs(matches_only = True, with_seq = True)
      second_pair_list = segment_reverse.get_aligned_pairs(matches_only = True, with_seq = True)
      pair_list = first_pair_list + second_pair_list
      chrom_start_position_int = pair_list[0][1] # 作为基准

      for tu in pair_list: # (0, 46232321, 'G')

         cycle = tu[1] - chrom_start_position_int

         if tu[2] in 'ATCG':
            matched_list.append(cycle)
         elif tu[2] in 'atcg':
            mismatch_list.append(cycle)
         else:
            continue

   elif segment_forward.flag == 163 and segment_reverse.flag == 83:  #  F2R1

      first_pair_list = segment_reverse.get_aligned_pairs(matches_only = True, with_seq = True)
      second_pair_list = segment_forward.get_aligned_pairs(matches_only = True, with_seq = True)
      pair_list = first_pair_list[::-1] + second_pair_list[::-1]
      chrom_start_position_int = pair_list[0][1] # 作为基准

      for tu in pair_list:
         cycle = chrom_start_position_int - tu[1]

         if tu[2] in 'ATCG':
            matched_list.append(cycle)
         elif tu[2] in 'atcg':
            mismatch_list.append(cycle)
         else:
            continue


   return matched_list, mismatch_list



def main(argvList = sys.argv, argv_int = len(sys.argv)):

   sam_file_str = path.realpath(path.expanduser(argvList[1]))
   output_str = 'summary_{}.txt'.format(path.splitext(path.basename(sam_file_str))[0])  # 如果sam_file_str为'/home/user/server/personal/test.sam'， output_str为summary_test.txt

   count_dict: dict[int:list] = {}
   read_pair_count_int = 0
   for forward_read, reverse_read in iter_queryname_sam(sam_file_str):
      match_list, mismatch_list = get_base_PCR_cycle(forward_read, reverse_read)
      for cycle_int in match_list:
         if cycle_int not in count_dict.keys():
            count_dict[cycle_int] = [0, 0]

         count_dict[cycle_int][1] += 1

      for cycle_int in mismatch_list:
         if cycle_int not in count_dict.keys():
            count_dict[cycle_int] = [0, 0]

         count_dict[cycle_int][0] += 1
         count_dict[cycle_int][1] += 1

      read_pair_count_int += 1
      if read_pair_count_int % 100000 == 0:
         keys_list = list(count_dict.keys())
         keys_list.sort()
         out_f = open(output_str, 'wt')
         for cycle_int in keys_list:
            line_str = f"Cycle: {cycle_int} {count_dict[cycle_int][0]}/{count_dict[cycle_int][1]} ({round(count_dict[cycle_int][0] / count_dict[cycle_int][1] * 1000, 4)}‰)"
            out_f.writelines(line_str + '\n')
         out_f.close()
         print(path.splitext(path.basename(sam_file_str))[0], read_pair_count_int, flush = True)

   keys_list = list(count_dict.keys())
   keys_list.sort()
   out_f = open(output_str, 'wt')
   for cycle_int in keys_list:
      line_str = f"Cycle: {cycle_int} {count_dict[cycle_int][0]}/{count_dict[cycle_int][1]} ({round(count_dict[cycle_int][0] / count_dict[cycle_int][1] * 1000, 4)}‰)"
      out_f.writelines(line_str + '\n')
   out_f.close()
   print(path.splitext(path.basename(sam_file_str))[0], read_pair_count_int, flush = True)
   print('Done', output_str)

   return

if __name__ == '__main__':
   main()
