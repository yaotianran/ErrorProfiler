#!/usr/bin/env python3
# v 1.3d
# ~/projects/noise/noise.py -b raw_data/nawgs0628B_sub360m_filter.bam raw_data/nawgs0628B_sub360m.queryname.sam

import sys
import os
import os.path as path
import argparse
import multiprocessing
import pysam
import collections

from multiprocessing import managers
import time
import queue

import lib.overlap as overlap
import lib.utils as utils

# INPUT, OUTPUT, MODE, COUNT, MAPQ, PROCESSES， WRITE_DICT_PER_OVERLAP
ARGUMENTS_DICT = {}

def get_arguments() -> int:
   '''
   获取命令行参数，返回一个字典

   Parameters:

   Returns:
       **arguments_dict**: dict
            命令行参数字典
   '''

   global ARGUMENTS_DICT

   argvList = sys.argv

   parser_ar = argparse.ArgumentParser(prog='PROG',
                                        usage=' python3 {0} [OPTION] INPUT'.format(argvList[0]),
                                        description='提取sam文件中pair-end的overlap区域的位点信息',
                                        epilog='epilog text',
                                        formatter_class=argparse.RawTextHelpFormatter)

   parser_ar.add_argument('INPUT', help= 'FILE. Queryname-sorted sam', metavar = 'Queryname-sorted sam file')

   parser_ar.add_argument('-o', '--output', default = '', help= '输出文件的前缀', metavar = '<result>', dest='OUTPUT')
   parser_ar.add_argument('-m', '--mode', default = 'raw', help= '模式～raw, bin. 默认为raw模式', metavar = '<raw>', dest='MODE')
   parser_ar.add_argument('-b', '--bam', default = '', help= 'sorted and index bam file. 如果提供了bam文件，则自动使用pysam判断真实碱基', metavar = '', dest='BAM')
   parser_ar.add_argument('-c', '--count', default = 1_000_000, type= int, help= 'INT. 检测overlap区域的数量', metavar = '<1_000_000>', dest= 'COUNT')
   parser_ar.add_argument('-q', '--min_mapq', default = 0, type= int, help= 'INT. 只检测不低于特定MAPQ的pair', metavar = '<0>', dest= 'MAPQ')
   parser_ar.add_argument('-l', '--max_tlen', default = 999, type= int, help= 'INT. 只检测不大于特定TLEN的pair (建议设置成读长+1)', metavar = '<999>', dest= 'TLEN')
   parser_ar.add_argument('-t', '--threads', default = 24, type= int, help= 'INT. 进程数', metavar = '<24>', dest= 'PROCESSES')


   paramters = parser_ar.parse_args()

   ARGUMENTS_DICT['INPUT'] = paramters.INPUT  # SAM file

   ARGUMENTS_DICT['OUTPUT'] = paramters.OUTPUT  #
   ARGUMENTS_DICT['MODE'] = paramters.MODE  #
   ARGUMENTS_DICT['BAM'] = paramters.BAM  #
   ARGUMENTS_DICT['COUNT'] = paramters.COUNT  #
   ARGUMENTS_DICT['MAPQ'] = paramters.MAPQ  #
   ARGUMENTS_DICT['TLEN'] = paramters.TLEN  #
   ARGUMENTS_DICT['PROCESSES'] = paramters.PROCESSES  #

   return 0

def __write_count_dict(count_dict, output_file, mode, minimum_mismatch: int = -1, minimum_total: int = -1, other = '', write_header = True) -> int:
   '''
   将count_dict写入指定文件

   Parameters:
      **count_dict**: {overlap.RawBin:int} 或者是 {overlap.Bin:[mismatch_base_count:int, total]}
         记录不同bin下面的base数量

      **output_file**: str
         写入的文件

   Returns:
       **value**: 0 on success

   '''
   mismatch_total_int = 0
   based_total_int = 0
   big_bin_int = 0
   small_bin_int = 0
   other_str = other
   t = time.time()


   with open(output_file, 'at') as out_f:
      total_bin_count_int = len(count_dict)
      i = 0

      if mode == 'bin':
         if write_header:
            header_str = '#bin\ncontext\tcycle\tqual_upstream\tqual\tmismatch\ttotal\n'
            out_f.writelines(header_str)
         for assigned_bin, count_lst in count_dict.items():
            line_str = '{context}\t{cycle}\t{qual_upstream}\t{qual}\t{mismatch}\t{total}\n'.format(context = assigned_bin.context,
                                                                                                   cycle = assigned_bin.cycle,
                                                                                                   qual_upstream = assigned_bin.qual_upstream,
                                                                                                   qual = assigned_bin.qual,
                                                                                                   mismatch = count_lst[0],
                                                                                                   total = count_lst[1],
                                                                                                   )
            out_f.writelines(line_str)
            i += 1
            if i % 100_000 == 0:
               message = 'Writing bins ({}/{}, {}%)                 '.format(utils.count_to_unit(i),
                                                           utils.count_to_unit(total_bin_count_int),
                                                           round(i / total_bin_count_int * 100, 1)
                                                           )
               print(message, end = '\r')




      if mode == 'raw':
         if write_header:
            header_str = '#raw bin\nread_order\tcycle\tgc_ratio\tbase_real\tqual_upstream_2\tqual_upstream\tqual\tqual_downstream\tqual_downstream_2\tbase_upstream_2\tbase_upstream\tbase_called\tbase_downstream\tbase_downstream_2\tcount\n'
            out_f.writelines(header_str)
         for assigned_raw_bin, count_int in count_dict.items():
            line_str = '{read_order}\t{cycle}\t{gc_ratio}\t{base_real}\t{qual_upstream_2}\t{qual_upstream}\t{qual}\t{qual_downstream}\t{qual_downstream_2}\t{base_upstream_2}\t{base_upstream}\t{base_called}\t{base_downstream}\t{base_downstream_2}\t{count}\n'.format(read_order = assigned_raw_bin.read_order,
                                                                                                                                                                                                                                                                         cycle = assigned_raw_bin.cycle,
                                                                                                                                                                                                                                                                         gc_ratio = assigned_raw_bin.gc_ratio,
                                                                                                                                                                                                                                                                         base_real = assigned_raw_bin.base_real,

                                                                                                                                                                                                                                                                         base_upstream_2 = assigned_raw_bin.base_upstream_2,
                                                                                                                                                                                                                                                                         base_upstream = assigned_raw_bin.base_upstream,
                                                                                                                                                                                                                                                                         base_called = assigned_raw_bin.base_called,
                                                                                                                                                                                                                                                                         base_downstream = assigned_raw_bin.base_downstream,
                                                                                                                                                                                                                                                                         base_downstream_2 = assigned_raw_bin.base_downstream_2,

                                                                                                                                                                                                                                                                         qual_upstream_2 = assigned_raw_bin.qual_upstream_2,
                                                                                                                                                                                                                                                                         qual_upstream = assigned_raw_bin.qual_upstream,
                                                                                                                                                                                                                                                                         qual = assigned_raw_bin.qual,
                                                                                                                                                                                                                                                                         qual_downstream = assigned_raw_bin.qual_downstream,
                                                                                                                                                                                                                                                                         qual_downstream_2 = assigned_raw_bin.qual_downstream_2,

                                                                                                                                                                                                                                                                         count = count_int,
                                                                                                                                                                                                                                                                         )
            out_f.writelines(line_str)
            i += 1
            if i % 100_000 == 0:
               message = 'Writing bins ({}/{}, {}%)                     '.format(utils.count_to_unit(i),
                                                           utils.count_to_unit(total_bin_count_int),
                                                           round(i / total_bin_count_int * 100, 1)
                                                           )
               print(message, end = '\r')


   message = '[{}] Written {} bins in {} seconds'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                                                            total_bin_count_int,
                                                            round(time.time() - t, 1)
                                                            )
   print(message)
   return 0

def assign_bin_and_record(write_queue: managers.queue, output_file: str, mode: str):
   """
   从write_queue读入base_info list, 按照mode, 分组为Bin或者RawBin

   写入到指定位置

   Parameters:
      **write_queue**: type
         its function

      **mode**:str
         模式，功能如下：
         raw ~ 在write_queue写入根据BaseInfo判断好的RawBin类
         bin ~ 在write_queue直接写入BaseInfo，让write_queue的接收端统计接受到的BaseInfo

   Returns:
       **value**: None


   """
   SHOW_PROGRESS_PER_OVERLAP = 1_000
   WRITE_PER_BINS = 30_000_000_000_000_000
   HIGH_QUALITY = 26

   if mode == 'bin':
      count_dict = collections.defaultdict(lambda :[0, 0])  # {Bin: List[mismatch:int, total:int]}

   if mode == 'raw':
      count_dict = collections.defaultdict(int) # {RawBin: int}


   total_baseinfo_int = 0 # 记录接收的总baseinfo list数量，有的baseinfo list为空[]
   effect_baseinfo_int = 0 # 记录不为空的baseinfo list数量

   total_base_int = 0 # 记录总位点个数
   total_mismatch_base_int = 0 # 记录总错误位点个数
   total_mismatch_high_qual_base_int = 0 # 记录错误位点中高质量错误位点的个数

   is_write_header_bool = True # 写入表格文件时受否有header

   t = time.time()
   print('Ready to receive bases ...', flush = True)

   real_base_source_lst = [0, 0, 0, 0, 0, 0, 0]
   while True:
      try:
         base_info_lst = write_queue.get() #  一个list,每个元素是BaseInfo TypedDict

         if base_info_lst == '#done#':
            raise queue.Empty

         total_baseinfo_int += 1

         if base_info_lst != []:
            effect_baseinfo_int += 1
            total_base_int += len(base_info_lst)
            for base_info in base_info_lst:
               real_base_source_lst[-base_info['base_real_source']] += 1

               if base_info['base_called'] != 'N' and base_info['base_real'] != 'N' and base_info['base_called'] != base_info['base_real']:
                  total_mismatch_base_int += 1
                  if base_info['base_qual'] >= HIGH_QUALITY:
                     total_mismatch_high_qual_base_int += 1

               if mode == 'bin':
                  assigned_bin = overlap.assign_bin(base_info)
                  if assigned_bin is not None:
                     count_dict[assigned_bin][1] += 1
                     if base_info['base_called'] != base_info['base_real']:
                        count_dict[assigned_bin][0] += 1

               if mode == 'raw':
                  # update student set score=score+1 where id = 1
                  assigned_raw_bin = overlap.assign_raw_bin(base_info)
                  if assigned_raw_bin is not None:
                     count_dict[assigned_raw_bin] += 1


         if total_baseinfo_int % SHOW_PROGRESS_PER_OVERLAP == 0:
            time_float = time.time() - t
            t = time.time()
            if len(count_dict) % WRITE_PER_BINS != 0:
               bin_left_int = WRITE_PER_BINS - len(count_dict) % WRITE_PER_BINS
            else:
               bin_left_int = 0
            message = '[{t}] Bases:{base} Bins:{bin_count} Overlaps:{eff_overlap}/{total_overlap}({percent}%) Bpo:{bpo} Bpb:{bpb} Mismatch:{mismatch}({rate1}‰) Highqual:{highqual}({rate2}‰) Speed:{speed}'.format(t = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                                                                                                                                                                                                                                    bin_count = utils.count_to_unit(len(count_dict)),
                                                                                                                                                                                                                                    eff_overlap = utils.count_to_unit(effect_baseinfo_int),
                                                                                                                                                                                                                                    total_overlap = utils.count_to_unit(total_baseinfo_int),
                                                                                                                                                                                                                                    percent = round(effect_baseinfo_int / total_baseinfo_int * 100, 1),
                                                                                                                                                                                                                                    #bin_left = utils.count_to_unit(bin_left_int),
                                                                                                                                                                                                                                    base = utils.count_to_unit(total_base_int),
                                                                                                                                                                                                                                    bpb = round(total_base_int / len(count_dict), 1),
                                                                                                                                                                                                                                    bpo = round(total_base_int / total_baseinfo_int / 2, 1),
                                                                                                                                                                                                                                    mismatch = utils.count_to_unit(total_mismatch_base_int),
                                                                                                                                                                                                                                    rate1 = round(total_mismatch_base_int / total_base_int * 1000, 2),
                                                                                                                                                                                                                                    highqual = utils.count_to_unit(total_mismatch_high_qual_base_int),
                                                                                                                                                                                                                                    rate2 = round(total_mismatch_high_qual_base_int / total_base_int * 1000, 2),
                                                                                                                                                                                                                                    speed = round(SHOW_PROGRESS_PER_OVERLAP / time_float)
                                                                                                                                                                                                                                    )
            print(message, flush = True)
            print(real_base_source_lst)


         if len(count_dict) >= WRITE_PER_BINS or os.access('write', os.R_OK):
            os.rename('write', 'writing')
            r = __write_count_dict(count_dict, output_file, mode, write_header = is_write_header_bool)
            is_write_header_bool = False
            del count_dict
            os.rename('writing', 'done_writen')


            if mode == 'bin':
               count_dict = collections.defaultdict(lambda :[0, 0])  # {Bin: List[mismatch:int, total:int]}

            if mode == 'raw':
               count_dict = collections.defaultdict(int) # {RawBin: int}


      except queue.Empty:
         time_float = time.time() - t
         message = '[{t}] Bases:{base} Bins:{bin_count} Overlaps:{eff_overlap}/{total_overlap}({percent}%) Bpo:{bpo} Bpb:{bpb} Mismatch:{mismatch}({rate1}‰) Highqual:{highqual}({rate2}‰) Speed:{speed}'.format(t = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                                                                                                                                                                                                                                 bin_count = utils.count_to_unit(len(count_dict)),
                                                                                                                                                                                                                                 eff_overlap = utils.count_to_unit(effect_baseinfo_int),
                                                                                                                                                                                                                                 total_overlap = utils.count_to_unit(total_baseinfo_int),
                                                                                                                                                                                                                                 percent = round(effect_baseinfo_int / total_baseinfo_int * 100, 1),
                                                                                                                                                                                                                                 #bin_left = utils.count_to_unit(bin_left_int),
                                                                                                                                                                                                                                 base = utils.count_to_unit(total_base_int),
                                                                                                                                                                                                                                 bpb = round(total_base_int / len(count_dict), 1),
                                                                                                                                                                                                                                 bpo = round(total_base_int / total_baseinfo_int / 2, 1),
                                                                                                                                                                                                                                 mismatch = utils.count_to_unit(total_mismatch_base_int),
                                                                                                                                                                                                                                 rate1 = round(total_mismatch_base_int / total_base_int * 1000, 2),
                                                                                                                                                                                                                                 highqual = utils.count_to_unit(total_mismatch_high_qual_base_int),
                                                                                                                                                                                                                                 rate2 = round(total_mismatch_high_qual_base_int / total_base_int * 1000, 2),
                                                                                                                                                                                                                                 speed = round(SHOW_PROGRESS_PER_OVERLAP / time_float)
                                                                                                                                                                                                                                 )
         print(message, flush = True)
         print(real_base_source_lst)


         r = __write_count_dict(count_dict, output_file, mode, write_header = is_write_header_bool)
         print('Write queue empty', flush = True)
         del count_dict
         break

      except Exception as ex:
         print('assign_bin_and_record:', ex, flush = True)

   return

def multiple_processes_worker(read_for_str: str, read_rev_str: str, bam_file: str, write_queue: managers.queue, base_counter: managers.Value):
   '''
   用于multiple processing，将read1和read2的分析结果，通过write_queue传输到记录进程，写入文件

   Parameters:
      **bam_file**: string
         一个sorted bam 文件或者为''. 为空字符时不使用pysam

   Returns:
       **value**: None

   '''

   if bam_file != '':
      bam_af = pysam.AlignmentFile(bam_file)
   else:
      bam_af = None

   base_info_lst = overlap.get_base_info(read_for_str, read_rev_str, bam_af)

   try:
      bam_af.close()
   except:
      pass

   # 传送一次的开销大，所以以一个read pair形成的base info list传输。
   # 此list在read pair没有overlap区域的情况下为[]，但是也可以传入write_queue，这样可以让assign_bin_and_record函数统计有效read pair的百分比
   write_queue.put(base_info_lst)

   #if len(base_info_lst) > 0:
      #base_counter.value += len(base_info_lst)

      # print(base_info_lst[0])
      # print(base_info_lst[1])
      # print('multiple_processes_worker send', len(base_info_lst), 'items')

   return


def main(argvList = sys.argv, argv_int = len(sys.argv)):

   # QUERYNAME_SAM
   QUERYNAME_SAM = path.realpath(path.expanduser(ARGUMENTS_DICT['INPUT']))
   if not os.access(QUERYNAME_SAM, os.R_OK):
      message = '无法读取输入的文件(输入为{})'.format(QUERYNAME_SAM)
      print(message)
      sys.exit(1)

   # output_str
   OUTPUT = ARGUMENTS_DICT['OUTPUT']
   if OUTPUT == '':
      prefix = path.splitext(path.basename(QUERYNAME_SAM))[0]
   else:
      prefix = path.realpath(path.expanduser(ARGUMENTS_DICT['OUTPUT']))

   MODE = ARGUMENTS_DICT['MODE']
   if MODE == 'bin':
      suffix = '.bin.txt'
   else:
      suffix = '.rawbin.txt'

   output_str = prefix + suffix
   if os.access(output_str, os.R_OK) and os.access(output_str, os.W_OK):
      os.remove(output_str)

   # SORTED_BAM
   BAM = ARGUMENTS_DICT['BAM']
   SORTED_BAM = path.realpath(path.expanduser(BAM))
   BAI = SORTED_BAM + '.bai'
   if BAM == '':
      SORTED_BAM = ''
   elif not os.access(SORTED_BAM, os.R_OK):
      message = '无法读取BAM文件(输入为{})'.format(SORTED_BAM)
      print(message)
      sys.exit(1)
   elif not os.access(BAI, os.R_OK):
      message = 'BAM文件没有index，使用pysam index...'
      print(message)
      r = utils.index(SORTED_BAM, PROCESSES)
   else:
      pass


   PROCESSES = ARGUMENTS_DICT['PROCESSES']
   MAPQ = ARGUMENTS_DICT['MAPQ']
   TLEN = ARGUMENTS_DICT['TLEN']
   COUNT = ARGUMENTS_DICT['COUNT']


   # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   t = time.time();
   manager = multiprocessing.Manager()
   base_counter = manager.Value('i', 0) # base_counter 用于记录收集的overlap区域base数量
   message = 'MODE:{} MIN_MAPQ:{} MAX_TLEN:{} OVERLAP_COUNT:{}'.format(MODE, MAPQ, TLEN, COUNT)
   print(message)

   # 首先启动“分类和写入”进程
   write_queue = manager.Queue()
   write_pool = multiprocessing.Pool(1)
   write_pool.apply_async(assign_bin_and_record, args = (write_queue, output_str, MODE, ))

   # use chunk pool method。
   # 常规pool占用大量内存。这里首先建立一个pool，使用该pool处理一个CHUNK，处理完当前CHUNK后会销毁这个pool释放内存，然后建立一个新pool。
   CHUNK_SIZE = PROCESSES * 10000      # 每个PROCESS处理10000个job，然后销毁
   jobs_lst = []
   i = 0  #  记录已处理的reads对的数量      # time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
   worker_pool = multiprocessing.Pool(PROCESSES)
   for read_for_str, read_rev_str in overlap.iter_queryname_sam(QUERYNAME_SAM, mapq = MAPQ, tlen = TLEN):
      if i >= COUNT:
         break

      job = worker_pool.apply_async(multiple_processes_worker, args = (read_for_str, read_rev_str, SORTED_BAM, write_queue, base_counter))
      jobs_lst.append(job)
      i += 1

      if len(jobs_lst) % CHUNK_SIZE == 0 and len(jobs_lst) != 0:
         message = '[{}] Start processing a chunk.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
         print(message)
         for job in jobs_lst:
            job.get()

         worker_pool.close()
         worker_pool.join()
         jobs_lst = []
         worker_pool = multiprocessing.Pool(PROCESSES)


   # 如果jobs_lst还有没有处理的job，处理完
   if len(jobs_lst) != 0:
      for job in jobs_lst:
         job.get()

   worker_pool.close()
   worker_pool.join()

   write_queue.put('#done#')
   write_pool.close()
   write_pool.join()
   message = 'Processed {} reads in {} seconds'.format(i, round(time.time() - t))
   print(message)

   return


if __name__ == '__main__':
   r = get_arguments()
   main()


