#!/usr/bin/env python3
# v1.0a
# 过滤重叠长度为特定值的read pairs
# 输出 filtered read count / total read count, mismatch base count / total base count
# filter_overlap_length.py queryname.sam 150
import os
import os.path as path
import sys
import re

SOFTCLIP_COMPILE = re.compile(r'\d+(?=S)')
SOFTCLIP_LEFT_COMPILE = re.compile(r'^\d+(?=S)')
SOFTCLIP_RIGHT_COMPILE = re.compile(r'\d+(?=S$)')
MATCH_COMPILE = re.compile(r'\d+(?=M)')
MD_COMPILE = re.compile(r'(?<=MD:Z:)\w+')
MD_SPLIT_COMPILE = re.compile(r'\d+|[A-Z]')
# 确定是否overlap和overlap区域的绝对/相对坐标
# read1，read2：sam文件的两个read条目
# 无overlap则返回None
# overlap_pos_left_int, overlap_pos_right_int, overlap_length_int, read_for_overlap_left, read_for_overlap_right, read_rev_overlap_left, read_rev_overlap_right
def get_overlap_coordinate(read_for_str: str, read_rev_str: str) -> tuple[int, int, int, int, int, int, int]:
   '''
   提取两个read的overlap坐标和序列等信息

   Parameters:
      **read_for_str，read_rev_str**: str
         read pair的read string

   Returns:
      如果没有overlap，则返回None

      如果有overlap，则返回以下参数

      **overlap_pos_left_int, overlap_pos_right_int, overlap_length_int**: int，overlap区域左侧/右侧在基因组的位置和overlap区域的长度

         保证overlap_pos_left_int一定小于等于overlap_pos_right_int

         保证 overlap_length_int = overlap_pos_right_int - overlap_pos_left_int + 1


      **read_for_overlap_left，read_for_overlap_right**: int，forward read的overlap区域的相对坐标


         保证 read_for_lst[9][read_for_overlap_left:read_for_overlap_right]是forward read overlap区域的序列

         保证 read_for_overlap_right - read_for_overlap_left = overlap_length_int


      **read_rev_overlap_left，read_rev_overlap_right**: int，reverse read的overlap区域的相对坐标

         保证 read_rev_lst[9][read_rev_overlap_left:read_rev_overlap_right]是reverse read overlap区域的序列

         保证 read_rev_overlap_right - read_rev_overlap_left = overlap_length_int

   '''

   # 分割read字符串。read_for_lst, read_rev_lst
   read_for_lst = read_for_str.strip().split()
   read_rev_lst = read_rev_str.strip().split()


   # 统计两个read的map区域的长度。read_for_mapped_length_int，read_rev_mapped_length_int
   # 统计两个read的map区域的绝对坐标。 read_for_mapped_pos_left_int，read_for_mapped_pos_right_int, read_rev_mapped_pos_left_int，read_rev_mapped_pos_right_int
   read_for_match_lst = MATCH_COMPILE.findall(read_for_lst[5])
   read_rev_match_lst = MATCH_COMPILE.findall(read_rev_lst[5])
   if read_for_match_lst == [] or read_rev_match_lst == []:  #  read没有map到基因组上
      return None

   read_for_mapped_length_int = int(read_for_match_lst[0])
   read_rev_mapped_length_int = int(read_rev_match_lst[0])

   read_for_mapped_pos_left_int = int(read_for_lst[3])
   read_for_mapped_pos_right_int = int(read_for_lst[3]) + read_for_mapped_length_int - 1

   read_rev_mapped_pos_left_int = int(read_rev_lst[3])
   read_rev_mapped_pos_right_int = int(read_rev_lst[3]) + read_rev_mapped_length_int -1

   # 判断是否overlap并计算overlap区域基因组坐标。overlap_pos_left_int，overlap_pos_right_int，overlap_length_int
   if read_for_mapped_pos_right_int < read_rev_mapped_pos_left_int or read_rev_mapped_pos_right_int < read_for_mapped_pos_left_int:     # 无overlap区域
      return None

   temp_lst = [read_for_mapped_pos_left_int, read_for_mapped_pos_right_int, read_rev_mapped_pos_left_int, read_rev_mapped_pos_right_int]
   temp_lst.sort()
   overlap_pos_left_int = temp_lst[1]  # overlap区域左侧基因组位置，保证为较小值
   overlap_pos_right_int = temp_lst[2]  # overlap区域右侧基因组位置，保证为较大值
   overlap_length_int = overlap_pos_right_int - overlap_pos_left_int + 1

   # 统计两个read的5'，3'端clip的长度。read_for_clip_left_length_int、read_for_clip_right_length_int、read_rev_clip_left_length_int、read_rev_clip_right_length_int
   clip_left_read_for_lst = SOFTCLIP_LEFT_COMPILE.findall(read_for_lst[5])
   clip_right_read_for_lst = SOFTCLIP_RIGHT_COMPILE.findall(read_for_lst[5])
   clip_left_read_rev_lst = SOFTCLIP_LEFT_COMPILE.findall(read_rev_lst[5])
   clip_right_read_rev_lst = SOFTCLIP_RIGHT_COMPILE.findall(read_rev_lst[5])

   if clip_left_read_for_lst == []:
      read_for_clip_left_length_int = 0
   else:
      read_for_clip_left_length_int = int(clip_left_read_for_lst[0])

   if clip_right_read_for_lst == []:
      read_for_clip_right_length_int = 0
   else:
      read_for_clip_right_length_int = int(clip_right_read_for_lst[0])

   if clip_left_read_rev_lst == []:
      read_rev_clip_left_length_int = 0
   else:
      read_rev_clip_left_length_int = int(clip_left_read_rev_lst[0])

   if clip_right_read_rev_lst == []:
      read_rev_clip_right_length_int = 0
   else:
      read_rev_clip_right_length_int = int(clip_right_read_rev_lst[0])


   # 计算overlap区域在两个read的内部的相对位置。read_for_overlap_left， read_for_overlap_right， read_rev_overlap_left， read_rev_overlap_right
   # 有两种overlap模式，取决于两个read的map区域的绝对坐标的排序(S表示softclip)：
   # 1. 正常模式，read_for_mapped_pos_left_int <= read_rev_mapped_pos_left_int <= read_for_mapped_pos_right_int <= read_rev_mapped_pos_right_int
   #  SSSS========================>>SSSSS
   #               SSSS<<====================SSSS

   # 2. 交错模式，这种情况比较少见。 read_rev_mapped_pos_left_int <= read_for_mapped_pos_left_int <= read_rev_mapped_pos_right_int <= read_for_mapped_pos_right_int
   #               SSS======================>>SSSS
   #  SSSS<<========================SSSS

   if read_rev_mapped_pos_left_int > read_for_mapped_pos_left_int: # 正常模式
      read_for_overlap_left = read_for_clip_left_length_int + (read_rev_mapped_pos_left_int - read_for_mapped_pos_left_int)
      read_for_overlap_right = read_for_overlap_left + overlap_length_int
      read_rev_overlap_left = read_rev_clip_left_length_int
      read_rev_overlap_right = read_rev_clip_left_length_int + overlap_length_int
   else:   #  交错模式
      read_for_overlap_left = read_for_clip_left_length_int
      read_for_overlap_right = read_for_clip_left_length_int + overlap_length_int
      read_rev_overlap_left = read_rev_clip_left_length_int + (read_for_mapped_pos_left_int - read_rev_mapped_pos_left_int)
      read_rev_overlap_right = read_rev_overlap_left + overlap_length_int

   return overlap_pos_left_int, overlap_pos_right_int, overlap_length_int, read_for_overlap_left, read_for_overlap_right, read_rev_overlap_left, read_rev_overlap_right


# 遍历一个queryname sorted sam文件，返回一对read string
# for forward_read, reverse_read in iter_queryname_sam(sam_file):
#     do samething about forward_read and reverse_read
#     line_str, 'header'
def iter_queryname_sam(sam_file: str, mapq = 0, tlen = 9999) -> tuple[str, str]:
   '''

   mapq: 只返回MAPQ值大于等于mapq的read pair

   确保输出read_for_str， read_rev_str

   '''

   samfile_str = path.realpath(path.expanduser(sam_file))
   with open(samfile_str) as in_f:

      for line_str in in_f:

         if line_str.startswith('@'):
            yield line_str, 'header'
            continue

         try:
            read1 = line_str
            read2 = next(in_f)

            read1_lst = read1.split('\t')
            read2_lst = read2.split('\t')

            while read1_lst[0] != read2_lst[0] or read1_lst[2] != read2_lst[2]: # queryname不同或者染色体不同，下移一个read
               read1 = read2
               read2 = next(in_f)
               read1_lst = read1.split('\t')
               read2_lst = read2.split('\t')

            if int(read1_lst[4]) < mapq or int(read1_lst[4]) < mapq: # MAPQ 太低
               continue

            if 'D' in read1_lst[5] or 'I' in read1_lst[5] or 'D' in read2_lst[5] or 'I' in read2_lst[5]: # cigar内有插入缺失
               continue

            if abs(int(read1_lst[8])) > tlen or abs(int(read2_lst[8])) > tlen: # TLEN  模板长度太长
               continue

            if read1_lst[1] == '99' and read2_lst[1] == '147' or read1_lst[1] == '163' and read2_lst[1] == '83':  # read1为forward，read2为reverse
               yield read1, read2
               continue

            if read1_lst[1] == '147' and read2_lst[1] == '99' or read1_lst[1] == '83' and read2_lst[1] == '163': # read1为reverse，read2为forward
               yield read2, read1
               continue
         except StopIteration as ex:
            print('iter_queryname_sam: 已读取全部文件')
            return


         except Exception as ex:
            print('iter_queryname_sam:', ex, '\n', read1, '\n', read2)
            continue

   return

def main(argvList = sys.argv, argv_int = len(sys.argv)):

   input_str = path.realpath(path.expanduser(argvList[1]))
   overlap_cutoff_int = int(argvList[2])


   output_str = path.basename(path.splitext(input_str)[0]) + f'_{overlap_cutoff_int}' + path.splitext(input_str)[1]



   m = 0
   n = 0
   total_mismatch_int = 0
   total_readLen_int = 0
   with open(output_str, 'wt', buffering = 2) as out_f:
      for forward_read, reverse_read in iter_queryname_sam(input_str):
         if reverse_read == 'header':
            out_f.writelines(forward_read)
            continue

         r = get_overlap_coordinate(forward_read, reverse_read)
         if r is None:  # no overlap
            overlap_length_int = 0
         else:  # 有overlap
            overlap_pos_left_int, overlap_pos_right_int, overlap_length_int, read_for_overlap_left, read_for_overlap_right, read_rev_overlap_left, read_rev_overlap_right = r

         if overlap_length_int >= overlap_cutoff_int:
            out_f.writelines(forward_read)
            out_f.writelines(reverse_read)
            n += 1

            MD_str = MD_COMPILE.findall(forward_read)[0]  #  '117G28'
            mismatch_int = len(re.findall(r'[ATCG]', MD_str))
            readLen_int = int(re.findall(r'\d+(?=M)', forward_read.split()[5])[0])
            total_mismatch_int += mismatch_int
            total_readLen_int += readLen_int

            MD_str = MD_COMPILE.findall(reverse_read)[0]  #  '117G28'
            mismatch_int = len(re.findall(r'[ATCG]', MD_str))
            readLen_int = int(re.findall(r'\d+(?=M)', reverse_read.split()[5])[0])
            total_mismatch_int += mismatch_int
            total_readLen_int += readLen_int

         m += 1
         if m % 100000 == 0:
            print(f"{n}/{m} ({round(n / m * 100, 1)}%)   {total_mismatch_int}/{total_readLen_int} ({round(total_mismatch_int / total_readLen_int * 1000, 3)}‰)        ", end = '\r', flush = True)


   print(f"{n}/{m} ({round(n / m * 100, 1)}%)   {total_mismatch_int}/{total_readLen_int} ({round(total_mismatch_int / total_readLen_int * 1000, 3)}‰)        ", flush = True)
   print('Done', output_str, flush = True)

   return


if __name__ == '__main__':
   print()
   main()
