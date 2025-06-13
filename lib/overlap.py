# 处理overlap函数
import os.path as path
from typing import TypedDict
import collections
import re
import pysam
from multiprocessing import managers
import multiprocessing
from time import time
import sys

SOFTCLIP_COMPILE = re.compile(r'\d+(?=S)')
SOFTCLIP_LEFT_COMPILE = re.compile(r'^\d+(?=S)')
SOFTCLIP_RIGHT_COMPILE = re.compile(r'\d+(?=S$)')
MATCH_COMPILE = re.compile(r'\d+(?=M)')
MD_COMPILE = re.compile(r'(?<=MD:Z:)\w+')
MD_SPLIT_COMPILE = re.compile(r'\d+|[A-Z]')


# 只用于记录base和质控，不用于给质量建模
RawBin = collections.namedtuple('RawBin', ['read_order',  #  'f' or 's'
                                           'cycle',
                                           'gc_ratio',
                                           'base_real',

                                           'qual_upstream_2', # 如果无上游2bp，则设置成-1
                                           'qual_upstream',   # 如果无上游1bp，则设置成-1
                                           'qual',
                                           'qual_downstream',   # 如果无下游1bp，则设置成-1
                                           'qual_downstream_2', # 如果无下游2bp，则设置成-1

                                           'base_upstream_2',   # 如果无上游2bp，则设置成''
                                           'base_upstream',   # 如果无上游1bp，则设置成''
                                           'base_called',
                                           'base_downstream',   # 如果无下游1bp，则设置成''
                                           'base_downstream_2', # 如果无下游2bp，则设置成''
                                           ],
                                defaults = [None] * 14)


# 所有base归类的bin，用于质量建模和Q值的校正
# cycle: int = None  # cycle ⊂ [1~150]。共150个
# context: str = None # 由called_base及其上游1bp所组成的2连碱基, 如果base是read的第一个碱基，则上游为‘’，如果该base在first read上则在前面加‘f’前缀， 如果是second read上则在前面加‘s’前缀
# context ∈ {fAA, fAT, fAC, fAG, ......,fA,fT,fC,fG, sAA, sAT, sAC, sAG, ......,sA,sT,sC,sG}
# 所以context共有(16 + 4)*2 = 40个
# qual: int # qual ⊂ [0, 33]。 共34个
# qual_upstream：int # qual ⊂ [-1, 0, 33]。（如果base是read的第一个碱基，则为-1）
# 理论上共有 150 x 40 x 34 x 34 = 6,936,000 bins
Bin = collections.namedtuple('Bin', ['cycle',
                                     'context',
                                     'qual_upstream',
                                     'qual', ],
                             defaults = [None] * 4)


class BaseInfo(TypedDict):
   base_real: str = None # called base所在位置的真实碱基种类，base_real ∈ {A，T，C，G, ''， N} # '' ~ deletion，N ~ 未知
   base_real_source: int = None # 0, -1, -2, -3, -4, -5, -6 用于记录real base的计算来源，base_real_source ∈ {-2, -1, >=0}。-2为从BAM的majority allele得到，-1为reference genome base，0为read1和read2相同，>0为read1和read2不同时，对方base的quality

   read_flag: str = None # read_flag ∈ {'83'，'163', '99', '147'}
   gc_ratio: int = None # gc_ratio ⊂ [0, 1]
   cycle: int = None  # cycle ⊂ [1, 150]
   base_called: str = None # called base的碱基种类，base_called ∈ {A，T，C，G, N}
   base_upstream: str = None #      记录上游碱基。如果base是第一个碱基（cycle为1），则为''
   base_downstream: str = None #    记录下游碱基。如果base是最后一个碱基（cycle为最大），则为''
   base_qual: int = None # base_qual ⊂ [0, 49]
   base_qual_upstream: list[int] = None #  记录上游碱基质量。base_qual_upstream ⊂ [-1, 49]  如果没有上游序列,则设置成[]
   base_qual_downstream: list[int] = None #  记录下游碱基质量。base_qual_downstream ⊂ [-1, 49], 如果没有下游序列,则设置成[]

   chrom: str = None #
   pos: int = None
   refer_base: str = None #  此处reference 的base 种类
   query_name = None

# 输入一个质量字符串,输出一个integer list, 如果输入为空字符,返回[]
def __str_to_qual(char: str, phred = 33) -> list[int, ...]:
   if char == '':
      return []

   return list(map(lambda x: ord(x) - phred, char))

# 输入一个DNA序列, 输出反向互补序列, 如果序列中有N,则原样输出
def __rev_comp(base: str) -> str:
   try:
      if len(base) == 0:
         return ''
      rev_comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', '': '', }
      return ''.join([rev_comp_dict[x] for x in base[::-1]])
   except:
      return None

# 从read string其中的MD tag和SEQ推断出map区域的reference sequence
# 输入一个read string，返回map区域的reference sequence
def __extract_reference_by_MD(read_str: str) -> str:
   '''
   read_str是一条read string， 内部必须有SEQ字段和MD tag字段

   返回这个read的map区域的reference sequence，无论read是正向还是反向read，reference sequence都以基因组正链方向显示

   如果有错误发生，则返回''
   '''
   if 'MD:Z' not in read_str:
      message = '__extract_reference_by_MD: read没有MG tag\n{}'.format(read_str)
      print(message)
      return ''

   read_lst = read_str.split('\t')

   clip_left_lst = SOFTCLIP_LEFT_COMPILE.findall(read_lst[5])  #  ['5']
   clip_left_length_int = int(clip_left_lst[0]) if clip_left_lst != [] else 0  # read 左侧soft clip的长度

   clip_right_lst = SOFTCLIP_RIGHT_COMPILE.findall(read_lst[5])  #  ['4']
   clip_right_length_int = int(clip_right_lst[0]) if clip_right_lst != [] else 0 # read 右侧soft clip的长度

   try:
      MD_str = MD_COMPILE.findall(read_str)[0]  #  '117G28'
      MD_lst = MD_SPLIT_COMPILE.findall(MD_str) #  ['117', 'G', '28']
      MD_length_int = sum([1 if x in 'ATCGN' else int(x) for x in MD_lst])

      mapped_end_left_int = len(read_lst[9]) - clip_right_length_int  #  mapped区域终止位置相对于序列左侧的坐标
      mapped_str = read_lst[9][clip_left_length_int:mapped_end_left_int]
      mapped_length_int = mapped_end_left_int - clip_left_length_int
      if mapped_length_int != MD_length_int:
         message = '__extract_reference_by_MD: {}(flag{}) mapped序列长度与MD tag长度不一致 (CIAGR:{} MD:{})'.format(read_lst[0], read_lst[1], read_lst[5], MD_str)
         print(message)
         raise IndexError
      else:
         reference_str = ''
         for s in MD_lst: #  ['117', 'G', '28']
            if s.isnumeric():
               reference_str += mapped_str[len(reference_str):(len(reference_str) + int(s))]
            elif s in 'ATCGN':
               reference_str += s
            else:
               message = '__extract_reference_by_MD: {}(flag{}) 非法MD (CIAGR:{} MD:{})'.format(read_lst[0], read_lst[1], read_lst[5], MD_str)
               print(message)
               raise IndexError

      if len(reference_str) != len(mapped_str):
         message = '__extract_reference_by_MD: {}(flag{}) mapped序列长度与MD tag长度不一致 (CIAGR:{} MD:{})'.format(read_lst[0], read_lst[1], read_lst[5], MD_str)
         print(message)
         raise IndexError

   except IndexError:
      reference_str = ''


   return reference_str


# 确定是否overlap和overlap区域的绝对/相对坐标
# read1，read2：sam文件的两个read条目
# 无overlap则返回None
# overlap_pos_left_int, overlap_pos_right_int, overlap_length_int, read_for_overlap_left, read_for_overlap_right, read_rev_overlap_left, read_rev_overlap_right
def __get_overlap_coordinate(read_for_str: str, read_rev_str: str) -> tuple[int, int, int, int, int, int, int]:
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



# 检索 bam 文件，得到chrom：pos位置的major allele
# 这里使用pysam，因为samtools mpileup命令运行更慢。这个函数须运行几亿次
# 返回ATCG之中的一个, 也有可能返回 '', 表示有两个或两个以上的major allele
def __get_major_allele(bam_af: pysam.AlignmentFile, chrom: str, pos: int) -> str:
   '''
   处理 bam 文件，得到相关位置的major allele
   这里使用pysam，因为samtools mpileup命令运行更慢
   这个函数须运行几亿次，要快

   Parameters:
       **bam_af**: pysam.AlignmentFile
         pysam.AlignmentFile object

       **chrom**: type
         染色体

      **pos**: int
         位置，兼容IGV，1-based

   Returns:
       **value**: str
         major allele
   '''
   from statistics import mean
   from itertools import compress

   for pileupcolumn in bam_af.pileup(contig = chrom,
                                     start=pos - 1,
                                     stop = pos,
                                     truncate = True,
                                     max_depth = 999999999,
                                     stepper = 'nofilter',
                                     ignore_overlaps = False,
                                     ignore_orphans = False,
                                     min_base_quality = 0):

      query_seq_lst = [x.upper() for x in pileupcolumn.get_query_sequences()]
      if len(query_seq_lst) < 4:  #  深度
         return None


      counter_dict = collections.Counter(query_seq_lst)

      base_lst = counter_dict.most_common() # [('A', 2), ('B', 2)]
      major_allele = ''
      if len(base_lst) == 1:  # 只有一个allele
         major_allele = base_lst[0][0]

      elif base_lst[0][1] > base_lst[1][1]: # 只有一个major allele
         major_allele = base_lst[0][0]

      else:
         pass

      return major_allele



      # 速度太慢，取消后续运行
      # # # # # # # 下方代码永远不会运行 # # # # # # #
      # 有两个或两个以上major allele, 按照平均测序质量选取
      base_candidate_lst = []
      for tu in base_lst:
         if tu[1] == base_lst[0][1]:
            base_candidate_lst.append(tu[0])

      query_qual_lst = pileupcolumn.get_query_qualities()

      major_allele = base_candidate_lst[0]
      fil = map(lambda x: True if x == major_allele else False, query_seq_lst)
      qual_mean_float = mean(compress(query_qual_lst, fil))

      for base in base_candidate_lst[1:]:
         fil = map(lambda x: True if x == base else False, query_seq_lst)
         q = mean(compress(query_qual_lst, fil))
         if q > qual_mean_float:
            major_allele = base
            qual_mean_float = q

      return major_allele

   return None

# 提取一个base的信息，返回一个BaseInfo dict (不提取gc_ratio, 不提取base_real，不提取base_real_source, 不提取refer_base)
# read_lst是一个read list。loc_from_read_left_end是相对于read左侧（5end）的坐标，0-based
# 反向序列的所有cycle，base和quality计算都基于反向互补序列
def __base_info(read_lst: list, loc_from_read_left_end: int, read_clip_left_length_int: int) -> BaseInfo:
   '''
   提取一个base的信息，返回一个BaseInfo dict

   Parameters:
      **read_lst**: list
         一个read list, 是原始sam文件条目，不作任何修改，分割得到

      **loc_from_read_left_end**: int
         该base相对于read左侧（5‘end）的坐标，0-based

   Returns:
      **BaseInfo**: TypedDict
         一个BaseInfo dict
   '''

   read_query_name_str = read_lst[0]
   read_flag_str = read_lst[1]
   read_chrom_str = read_lst[2]
   read_seq_str = read_lst[9]
   read_qual_str = read_lst[10]
   read_qual_integer_lst = __str_to_qual(read_qual_str)


   # genomic position
   pos = loc_from_read_left_end - read_clip_left_length_int + int(read_lst[3])


   # other attributes
   if read_flag_str == '99' or  read_flag_str == '163':  #  forward read
      cycle = loc_from_read_left_end + 1

      base_upstream = read_seq_str[:loc_from_read_left_end]
      base_called = read_seq_str[loc_from_read_left_end]
      base_downstream = read_seq_str[loc_from_read_left_end + 1:]

      base_qual_upstream = read_qual_integer_lst[:loc_from_read_left_end]
      base_qual = read_qual_integer_lst[loc_from_read_left_end]
      base_qual_downstream = read_qual_integer_lst[loc_from_read_left_end + 1:]

   elif read_flag_str == '147' or  read_flag_str == '83':  #  reverse read

      # 因为无论对于正向还是反向序列，SAM中SEQ字段储存的是都是正向序列
      # 所以我们先将SEQ字段转化为反向互补序列，
      # 以下所有cycle，base和quality计算都基于反向互补序列
      read_seq_rc_str = __rev_comp(read_seq_str)
      read_qual_integer_rc_lst = read_qual_integer_lst[::-1]
      loc_from_read_left_end_rc = len(read_seq_rc_str) - loc_from_read_left_end - 1 # 在反向互补序列中，该base相对于左端的位置（0-based）。例如序列长度10, base左端位置为2, cycle为3, 那么在该序列的反向互补序列中，新的左端位置为10 - 2 - 1 = 7， cycle为7 + 1 = 8

      cycle = loc_from_read_left_end_rc + 1

      base_upstream = read_seq_rc_str[:loc_from_read_left_end_rc]
      base_called = read_seq_rc_str[loc_from_read_left_end_rc]
      base_downstream = read_seq_rc_str[loc_from_read_left_end_rc + 1:]

      base_qual_upstream = read_qual_integer_rc_lst[:loc_from_read_left_end_rc]
      base_qual = read_qual_integer_rc_lst[loc_from_read_left_end_rc]
      base_qual_downstream = read_qual_integer_rc_lst[loc_from_read_left_end_rc + 1:]

   else:
      return None


   base_info_dict = BaseInfo(read_flag = read_flag_str,
                             cycle = cycle,
                             base_called = base_called,
                             base_upstream = base_upstream,
                             base_downstream = base_downstream,
                             base_qual = base_qual,
                             base_qual_upstream = base_qual_upstream,
                             base_qual_downstream = base_qual_downstream,
                             chrom = read_chrom_str,
                             pos = pos,
                             query_name = read_query_name_str)

   return base_info_dict

# 根据base的位点信息，将它分配到一个Bin
# 输入一个BaseInfo TypedDict，返回一个Bin namedtuple
# 如果不属于任何一个bin，则返回None
def assign_bin(base_info: BaseInfo) -> Bin:

   if base_info['base_called'] == 'N' or base_info['base_real'] == 'N':
      return None


   # cycle
   cycle = base_info['cycle']


   # context
   try:
      if base_info['base_upstream'] == '':
         context = base_info['base_called']
      else:
         context =  base_info['base_upstream'][-1] + base_info['base_called']
   except Exception as ex:
      print('assign_bin:', ex)
      print(base_info)
      return None

   if base_info['read_flag'] == '99' or base_info['read_flag'] == '83':
      context = 'f' + context
   elif base_info['read_flag'] == '163' or base_info['read_flag'] == '147':
      context = 's' + context
   else:
      return None


   # qual_upstream
   try:
      qual_upstream = base_info['base_qual_upstream'][-1]
   except IndexError:
      qual_upstream = -1


   # qual
   qual = base_info['base_qual']

   assigned_bin = Bin(cycle = cycle,
                      context = context,
                      qual_upstream = qual_upstream,
                      qual = qual
                      )

   return assigned_bin

# 根据base的位点信息，将它分配到一个RawBin
# 输入一个BaseInfo TypedDict，返回一个RawBin namedtuple
# 如果不属于任何一个bin，则返回None
def assign_raw_bin(base_info: BaseInfo) -> RawBin:

   if base_info['base_called'] == 'N' or base_info['base_real'] == 'N':
      return None


   # = = = = = = = = read_flag = = = = = = = =
   if base_info['read_flag'] in ['99', '83']:
      read_order = 'f'
   elif base_info['read_flag'] in ['163', '147']:
      read_order = 's'
   else:
      read_order = ''

   # = = = = = = = = cycle = = = = = = = =
   cycle = base_info['cycle']


   # = = = = = = = = gc_ratio = = = = = = = =
   if base_info['gc_ratio'] <= 1 / 3:
      gc_ratio = 'low'
   elif 1 / 3 < base_info['gc_ratio'] and base_info['gc_ratio'] < 2 / 3:
      gc_ratio = 'medium'
   else:
      gc_ratio = 'high'

   # = = = = = = = = base read = = = = = = = =
   base_real = base_info['base_real']

   # = = = = = = = = qual = = = = = = = =
   try:
      qual_upstream_2 = base_info['base_qual_upstream'][-2]
   except IndexError:
      qual_upstream_2 = -1

   try:
      qual_upstream = base_info['base_qual_upstream'][-1]
   except IndexError:
      qual_upstream = -1

   qual = base_info['base_qual']

   try:
      qual_downstream = base_info['base_qual_downstream'][0]
   except IndexError:
      qual_downstream = -1

   try:
      qual_downstream_2 = base_info['base_qual_downstream'][1]
   except IndexError:
      qual_downstream_2 = -1

   # = = = = = = = = base = = = = = = = =
   try:
      base_upstream_2 = base_info['base_upstream'][-2]
   except IndexError:
      base_upstream_2 = ''

   try:
      base_upstream = base_info['base_upstream'][-1]
   except IndexError:
      base_upstream = ''

   base_called = base_info['base_called']

   try:
      base_downstream = base_info['base_downstream'][0]
   except IndexError:
      base_downstream = ''

   try:
      base_downstream_2 = base_info['base_downstream'][1]
   except IndexError:
      base_downstream_2 = ''

   # = = = = = = = = 赋值 = = = = = = = =
   assigned_raw_bin = RawBin(read_order = read_order,
                             cycle = cycle,
                             gc_ratio = gc_ratio,
                             base_real = base_real,

                             qual_upstream_2 = qual_upstream_2,
                             qual_upstream = qual_upstream,
                             qual = qual,
                             qual_downstream = qual_downstream,
                             qual_downstream_2 = qual_downstream_2,

                             base_upstream_2 = base_upstream_2,
                             base_upstream = base_upstream,
                             base_called = base_called,
                             base_downstream = base_downstream,
                             base_downstream_2 = base_downstream_2,
                             )

   return assigned_raw_bin


# 遍历一个queryname sorted sam文件，返回一对read string
# for read1, read2 in iter_queryname_sam(sam_file):
#     do samething about read1 and read2
def iter_queryname_sam(sam_file: str, mapq = 0, tlen = 300) -> tuple[str, str]:
   '''

   mapq: 只返回MAPQ值大于等于mapq的read pair

   确保输出read_for_str， read_rev_str

   '''

   samfile_str = path.realpath(path.expanduser(sam_file))
   with open(samfile_str) as in_f:

      for line_str in in_f:

         if line_str.startswith('@'):
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

# 输入一对read pair，获取overlap区域一系列位点信息
# 输入一个read1 , read2, bam_file。bam_file用来获取某一位置的major allele
# 返回一个list,每个元素是BaseInfo TypedDict
def get_base_info(read_for_str: str, read_rev_str: str, bam_af: pysam.AlignmentFile, q_cutoff = 25, q_diff: int = 10) -> list[BaseInfo, ...]:
   '''
   提取两个read的overlap中每一个base的信息，每个base的信息对应一个BaseInfo类实例
   提取gc_ratio, 提取base_real，提取base_real_source, 提取refer_base

   Parameters:
      **read1, read2**: str
         一对read string

      **bam_af**: pysam.AlignmentFile
         该样本的pysam.AlignmentFile，用来获取mismatch位点的major allele。
         如果设置为None，将不使用pysam获取major_allele。由于获取major allele是一个低速过程，在mismatch较多的样本中，获取major allele将降低运行速度

      **q_cutoff, q_diff**: int
         由于获取major allele是一个慢速过程，所以当mismatch位点中较高质量的Q值大于q_cutoff， 且两者差异大于等于q_diff时，我们直接使用Q值较大的位点作为真实位点。

      **mode**:str
         如果为‘all’，则输出所有位点信息
         如果为其他，则只输出mismatch位点信息

   Returns:
      **base_info_list**: list[BaseInfo, ...]
         每个元素是一个BaseInfo类实例，代表一个碱基的信息
   '''

   result_lst = []

   # 不使用__rev_comp以加快速度
   rev_comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', '': '', }

   read_for_lst = read_for_str.strip().split()
   read_rev_lst = read_rev_str.strip().split()

   # print(f'querying {read_for_lst[0]}')

   gc_for_float = sum([s in 'CG' for s in read_for_lst[9]]) / len(read_for_lst[9])
   gc_rev_float = sum([s in 'CG' for s in read_rev_lst[9]]) / len(read_rev_lst[9])

   clip_left_read_for_lst = SOFTCLIP_LEFT_COMPILE.findall(read_for_lst[5])
   clip_left_read_for_length_int = int(clip_left_read_for_lst[0]) if clip_left_read_for_lst != [] else 0


   clip_left_read_rev_lst = SOFTCLIP_LEFT_COMPILE.findall(read_rev_lst[5])
   clip_left_read_rev_length_int = int(clip_left_read_rev_lst[0]) if clip_left_read_rev_lst != [] else 0


   # = = = = = = = = = = = = 获取overlap区域坐标, 并遍历位点，获取信息 = = = = = = = = = = = =
   r = __get_overlap_coordinate(read_for_str, read_rev_str)
   if r is None:  # no overlap
      return []
   else:
      # 获取overlap区域
      overlap_pos_left_int, overlap_pos_right_int, overlap_length_int, read_for_overlap_left, read_for_overlap_right, read_rev_overlap_left, read_rev_overlap_right = r

      # 获取map区域的reference seqence
      reference_for_str = __extract_reference_by_MD(read_for_str)
      reference_rev_str = __extract_reference_by_MD(read_rev_str)

   # 不必对__get_overlap_coordinate的返回值做长度检查，返回值的正确性由函数__get_overlap_coordinate本身保证
   for i in range(overlap_length_int):

      base_pos_int = overlap_pos_left_int + i  # 该位点在基因组的位置
      base_left_read_for_int = read_for_overlap_left + i # 该位点相对于正向序列左侧位置, 0-based
      base_left_read_rev_int = read_rev_overlap_left + i # 该位点相对于反向序列左侧位置

      base_for_base_info = __base_info(read_for_lst, base_left_read_for_int, clip_left_read_for_length_int)
      base_rev_base_info = __base_info(read_rev_lst, base_left_read_rev_int, clip_left_read_rev_length_int)

      base_for_base_info['gc_ratio'] = gc_for_float
      base_rev_base_info['gc_ratio'] = gc_rev_float

      base_for_base_info['refer_base'] = reference_for_str[base_left_read_for_int - clip_left_read_for_length_int]
      base_rev_base_info['refer_base'] = rev_comp_dict[(reference_rev_str[base_left_read_rev_int - clip_left_read_rev_length_int])]

      if base_for_base_info['base_called'] == 'N' or base_rev_base_info['base_called'] == 'N' or base_for_base_info['refer_base'] == 'N' or base_rev_base_info['refer_base'] == 'N':
         continue


      # 判断真实碱基
      if base_for_base_info['base_called'] == rev_comp_dict[base_rev_base_info['base_called']]:   # 这是一个match碱基，called base即为真实碱基

         base_for_base_info['base_real'] = base_for_base_info['base_called']
         base_rev_base_info['base_real'] = base_rev_base_info['base_called']
         base_for_base_info['base_real_source'] = 0
         base_rev_base_info['base_real_source'] = 0


      elif base_for_base_info['base_called'] == base_for_base_info['refer_base'] :  # 这是一个mismatch，正向碱基与它的refer base相符，所以正向碱基为正确碱基

         base_for_base_info['base_real'] = base_for_base_info['base_called']
         base_rev_base_info['base_real'] = rev_comp_dict[base_for_base_info['base_called']]
         base_for_base_info['base_real_source'] = 0
         base_rev_base_info['base_real_source'] = -1

         # print(read_for_lst[0], base_for_base_info['pos'])

      elif base_rev_base_info['base_called'] == base_rev_base_info['refer_base']:  # 这是一个mismatch， 反向碱基与它的refer base相符，所以反向碱基为正确碱基

         base_for_base_info['base_real'] = rev_comp_dict[base_rev_base_info['base_called']]
         base_rev_base_info['base_real'] = base_rev_base_info['base_called']
         base_for_base_info['base_real_source'] = -2
         base_rev_base_info['base_real_source'] = 0

         # print(read_for_lst[0], base_for_base_info['pos'])

      elif base_for_base_info['base_qual'] >= q_cutoff and base_for_base_info['base_qual'] - base_rev_base_info['base_qual'] >= q_diff:  # 这是一个mismatch，正向碱基质量高，所以正向碱基为正确碱基

         base_for_base_info['base_real'] = base_for_base_info['base_called']
         base_rev_base_info['base_real'] = rev_comp_dict[base_for_base_info['base_called']]
         base_for_base_info['base_real_source'] = 0
         base_rev_base_info['base_real_source'] = -3

         # print(read_for_lst[0], base_for_base_info['pos'])

      elif base_rev_base_info['base_qual'] >= q_cutoff and base_rev_base_info['base_qual'] - base_for_base_info['base_qual'] >= q_diff: # 这是一个mismatch，反向碱基质量高，

         base_for_base_info['base_real'] = rev_comp_dict[base_rev_base_info['base_called']]
         base_rev_base_info['base_real'] = base_rev_base_info['base_called']
         base_for_base_info['base_real_source'] = -4
         base_rev_base_info['base_real_source'] = 0

         # print(read_for_lst[0], base_for_base_info['pos'])

      else:  # 使用pysam获取major allele
         try:
            # base_real = 'A'
            base_real = __get_major_allele(bam_af, base_for_base_info['chrom'], base_pos_int) # 可能为‘’，表示deletion

            if base_real in 'ATCG':
               base_for_base_info['base_real'] = base_real
               base_rev_base_info['base_real'] = rev_comp_dict[base_real]

               if base_for_base_info['base_called'] == base_real:
                  base_for_base_info['base_real_source'] = 0
               else:
                  base_for_base_info['base_real_source'] = -5

               if base_rev_base_info['base_called'] == rev_comp_dict[base_real]:
                  base_rev_base_info['base_real_source'] = 0
               else:
                  base_rev_base_info['base_real_source'] = -5


            else:   # 无法获取 real base
               base_for_base_info['base_real'] = 'U'
               base_rev_base_info['base_real'] = 'U'
               base_for_base_info['base_real_source'] = -6
               base_rev_base_info['base_real_source'] = -6

               message = '使用pysam获取到位置{}:{}的major allele为{}，不属于ATCG, 真实碱基将被标记成U'.format(base_for_base_info['chrom'], base_pos_int, base_real)
               print(message)


         except Exception as ex:
            base_for_base_info['base_real'] = 'U'
            base_rev_base_info['base_real'] = 'U'
            base_for_base_info['base_real_source'] = -6
            base_rev_base_info['base_real_source'] = -6

            message = 'get_base_info: {} {}:{} {}'.format(base_for_base_info['query_name'], base_for_base_info['chrom'], base_pos_int, ex)
            print(message)
            continue


      result_lst.append(base_for_base_info)
      result_lst.append(base_rev_base_info)


   return result_lst

# 读入一个report文件，tsv格式，with header，返回一个bin_dict。
# report文件必须有以下五列 ['context', 'qual', 'cycle', 'qual_upstream', 'log_error_rate_fit']
# 其中context 类型为str， context ∈ {fAA, fAT, fAC, fAG, ......,fA,fT,fC,fG, sAA, sAT, sAC, sAG, ......,sA,sT,sC,sG}， 共有(16 + 4)*2 = 40个
# qual 类型为 非负整数 # qual ⊂ [0, ....]。
# cycle为string，格式为 \d,\d...., 例如 1,2 或者134
# qual_upstream为string，格式为 \d,\d...., 例如 1,2 或者16
# log_error_rate_fit为浮点数

# 输出为 report_dict，其格式为{merged_Bin: str, ...}
# merged_Bin与bin格式基本相同，不同之处在于cycle和qual_upstream字段为[int, ...]
def read_report(report_file: str) -> dict[Bin: [str, float]]:

   header_lst = ['context', 'qual', 'cycle', 'qual_upstream', 'log_error_rate_fit']

   report_dict = {}

   report_file_str = path.realpath(path.expanduser(report_file))
   with open(report_file_str) as in_f:
      for i, line_str in enumerate(in_f):
         line_lst = line_str.strip().split()

         if i == 0:
            index_lst = list(map(lambda x: line_lst.index(x), header_lst)) #  must convert to list
            continue

         if i % 10000 == 0:
            print(i)

         try:
            temp_lst = [line_lst[x] for x in index_lst] #  temp_lst的顺序与header一致
         except ValueError as ex:   # 列不全
            message = 'read_report: report file 列不全 (第 {} 行：{})'.format(i, line_str)
            print(ex, '\n', message)
            sys.exit(-1)


         cycle_lst = [int(x) for x in temp_lst[2].split(',')]
         qual_upstream_lst = [int(x) for x in temp_lst[3].split(',')]

         try:
            bin_tu = Bin(context = temp_lst[0],
                         qual = int(temp_lst[1]),
                         cycle = cycle_lst,
                         qual_upstream = qual_upstream_lst
                        )

            log_error_rate_float = float(temp_lst[4])  # 防止 log_error_rate_float 溢出
            if log_error_rate_float > 0:
               log_error_rate_float = 0
            elif log_error_rate_float < -6:
               log_error_rate_float = -6
            else:
               pass

            letter_str = chr(round(-log_error_rate_float * 10 + 33))
            report_dict[bin_tu] = letter_str

         except ValueError as ex:  # 值不对
            message = 'read_report: report file 列不全 (第 {} 行：{})'.format(i, line_str)
            print(ex, '\n', message)
            sys.exit(-1)

   return report_dict

def __get_quality_letter(base_info: BaseInfo, report_dict: dict[Bin: str, ...]) -> str:


   # 发成错误返回空字符
   return



# 读入一组samfile的read pair，更改quality
# bin_qual_dict是read_report函数读入的bin dictionary
# read_for_str, read_rev_str = change_quality()
def change_quality(read1: str, read2: str, report_dict: dict[Bin: str, ...], quality_shift_for_unknown: int = 0, overlap_mismatch_correction: bool = False, major_allele_dict: dict = None) -> tuple[str, str]:

   '''
   进行两轮校正，第一轮分别按照bin_qual_dict，校正两个read的质量。第二轮遍历两个read的overlap区域，如果有不同的basecall，则按照major_allele_dict的更改，如果

   quality_shift_for_unknown 为当无法判断base属于哪个bin时，将原有质量字母，扩大或减小几个单位

   返回正向，反向read string

   '''

   read_for_str, read_rev_str = __check_pair(read1, read2)
   read_for_lst = read_for_str.strip().split()
   read_rev_lst = read_rev_str.strip().split()

   for_qual_ori_str = read_for_lst[10]
   rev_qual_ori_str = read_rev_lst[10]

   # 遍历forward read quality string
   for_qual_new_str = ''
   for i, original_quality_letter_char in enumerate(for_qual_ori_str):
      try:
         base_info = __base_info(read_for_lst, i)
         new_quality_letter_char = __get_quality_letter(base_info, report_dict)
         if new_quality_letter_char == '':
            new_quality_letter_char = chr(ord(original_quality_letter_char) + quality_shift_for_unknown)
      except Exception as ex:
         print('change_quality:', ex)
         new_quality_letter_char = original_quality_letter_char

      for_qual_new_str += new_quality_letter_char

   read_for_lst[10] = for_qual_new_str

   # 遍历reverse read  quality string
   rev_qual_new_str = ''
   for i, original_quality_letter_char in enumerate(rev_qual_ori_str):
      try:
         base_info = __base_info(read_rev_lst, i)
         new_quality_letter_char = __get_quality_letter(base_info, report_dict)
         if new_quality_letter_char == '':
            new_quality_letter_char = chr(ord(original_quality_letter_char) + quality_shift_for_unknown)

      except Exception as ex:
         print('change_quality:', ex)
         new_quality_letter_char = original_quality_letter_char

      rev_qual_new_str += new_quality_letter_char

   read_rev_lst[10] = rev_qual_new_str


   if overlap_mismatch_correction:
      # 获取overlap区域坐标
      r = __get_overlap_coordinate(read_for_str, read_rev_str)
      if r is None:  # no overlap, 直接返回
         read_for_str = '\t'.join(read_for_lst)
         read_rev_str = '\t'.join(read_rev_lst)
         return read_for_str, read_rev_str

      else:
         overlap_pos_left_int, overlap_pos_right_int, overlap_length_int, read_for_overlap_left, read_for_overlap_right, read_rev_overlap_left, read_rev_overlap_right = r

      # 第二轮遍历两个read的overlap区域，如果有不同的base call，则质量大者为real base call
      read_for_overlap_str = read_for_lst[9][read_for_overlap_left:read_for_overlap_right]
      read_rev_overlap_str = read_rev_lst[9][read_rev_overlap_left:read_rev_overlap_right]

      read_for_seq_lst = list(read_for_lst[9])
      read_rev_seq_lst = list(read_rev_lst[9])
      for i, tu in enumerate(zip(read_for_overlap_str, read_rev_overlap_str)):
         if tu[0] != tu[1]:
            base_left_read_for_int = read_for_overlap_left + i # 该位点相对于正向序列左侧位置, 0-based
            base_left_read_rev_int = read_rev_overlap_left + i # 该位点相对于反向序列左侧位置

            if for_qual_new_str[base_left_read_for_int] > rev_qual_new_str[base_left_read_rev_int]: #  forward read 质量高
               read_rev_seq_lst[base_left_read_rev_int] = read_for_seq_lst[base_left_read_for_int]
               # read_rev_qual_lst[base_left_read_rev_int] = read_for_qual_lst[base_left_read_for_int]

            elif rev_qual_new_str[base_left_read_rev_int] > for_qual_new_str[base_left_read_for_int]: #  reverse read 质量高
               read_for_seq_lst[base_left_read_for_int] = read_rev_seq_lst[base_left_read_rev_int]
               # read_for_qual_lst[base_left_read_for_int] = read_rev_qual_lst[base_left_read_rev_int]

            else:  # 质量相同
               base_info = __base_info(read_for_lst, base_left_read_for_int)
               assign_bin = assign_bin(base_info)
               if assign_bin in bin_qual_dict.keys():
                  read_for_qual_float = bin_qual_dict[assign_bin][1]
               else:
                  continue

               base_info = __base_info(read_rev_lst, base_left_read_rev_int)
               assign_bin = assign_bin(base_info)
               if assign_bin in bin_qual_dict.keys():
                  read_rev_qual_float = bin_qual_dict[assign_bin][1]
               else:
                  continue

               if read_for_qual_float > read_rev_qual_float:
                  read_rev_seq_lst[base_left_read_rev_int] = read_for_seq_lst[base_left_read_for_int]
                  # read_rev_qual_lst[base_left_read_rev_int] = read_for_qual_lst[base_left_read_for_int]

               elif read_rev_qual_float > read_for_qual_float:
                  read_for_seq_lst[base_left_read_for_int] = read_rev_seq_lst[base_left_read_rev_int]
                  # read_for_qual_lst[base_left_read_for_int] = read_rev_qual_lst[base_left_read_rev_int]

               else:
                  continue

      read_for_lst[9] = ''.join(read_for_seq_lst)
      read_rev_lst[9] = ''.join(read_rev_seq_lst)

   # read_for_lst[10] = ''.join(read_for_qual_lst)
   # read_rev_lst[10] = ''.join(read_rev_qual_lst)

   read_for_lst.append('XX:Z:{}'.format(for_qual_ori_str))
   read_rev_lst.append('XX:Z:{}'.format(rev_qual_ori_str))

   read_for_str = '\t'.join(read_for_lst)
   read_rev_str = '\t'.join(read_rev_lst)

   return read_for_str, read_rev_str

if __name__ == '__main__':
   i = 1
   '''
   # read1 = r"C2305160001:S:PRM32401030122:1:000000:R001:C019	83	12	133219225	60	146M4S	=	133219198	-173	GGAATTCCTCCAAGACAGGAATTTCACTGGCCAGCCTCTTCAGCTCCCAGCTGGACTGAACAGCGATGAGTGTGGGCCCCCGGCGCTCCTCCTGGGTCAAGGCAAAATGGAAGAAAANACCTGGGTGGACCCAGCCTCAAAGAAGACTTT	@A@AA@A2@B@A@@AA@AA?A@?<?BA@AAAAA?@A?A<@AA?AA@B0AAAAA@A@BA@@@AA@AAA=BA@A@AA=5@@@BA@@A@@AA@@A@@@@BBA@@A@BA@A?@B@A@@A@@!@A@@A@@@AAA@?B@@@AAB@ABB@AA@@;AA	NM:i:1	MD:Z:117G28	MC:Z:4S146M	AS:i:144	XS:i:0	RG:Z:sample2_salus_500M"
   # read2 = r"C2305160001:S:PRM32401030122:1:000000:R001:C019	163	12	133219198	60	4S146M	=	133219225	173	ACTGTGTCAGCCACACAGATAGGCACCAGTGGGAATTCCTCCAAGACAGGAATTTCACTGGCCAGCCTCTTCAGCTCCCAGCTGGACTGAACAGCGATGAGTGTGGGCCCCCGGCGCTCCTCCTGGGTCAAGGCAAAATGGAAGAAAAGA	@@@@@A@@@AAA@@B@@@A@?A@B@A@AA@@@@BAA@A@A@AAAAB@@@A@@BAA@@BBAA@@@AB@A@A@@?@@B@@@A@B@AA?@AA@@BAB@BA@@??A@@AA@AA@@AA@AB@?@A@@BB@@@@@ABA@@?BBA@@?AAA:A@A@@	NM:i:0	MD:Z:146	MC:Z:146M4S	AS:i:146	XS:i:0	RG:Z:sample2_salus_500M"
   # read1 = r"00000000000:S:00000000000000:1:000001:R004:C065	99	1	15932376	60	150M	=	15932400	174	GGCACAAAGTCAATATTTGGAAATGTTAGCAATCATTTGGAAGGTGCTTGACGCACAGCGTTGATCGATAAATGGGACCTATTCTTCTTGGTAGCACTGCTGGTCTAGAGCAGGATTCTGGTTTCTGAGCTCTCCCTGNNCNNAACNTCC	AA@BAB@@A@AAAABAAAA@BAABA>AAAAAAAA@AABBA@A@A?@A?BAAABA@AAABAAA@AB@A@@@@AAAAA@AABAAAAAA@AABABA?AAAAAAABAABBAA?A?ABA@BAA@A@A@@AAB@A?AB@A=<@A!!2!!B;B!=:?	NM:i:5	MD:Z:138C0T1T0C3C3	MC:Z:150M	AS:i:140	XS:i:19	RG:Z:na12878_wgs"
   # read = r"00000000000:S:00000000000000:1:000001:R004:C065	147	1	15932400	60	150M	=	15932376	-174	NNNNGCANTCATTTGGAAGGTGCTTGACGCACAGCGTTGATCGATAAATGGGACCTATTCTTCTTGGTAGCACTGCTGGTCTAGAGCAGGATTCTGGTTTCTGAGCTCTCCCTGCTCTCAACCTCCACTCCACCCCATCCTGACTTAAAA	!!!!:A@!@AA:A?>?@AB@?=@/@@BA>>AA=@BAA?A??@?BA7<?AAB??AB?AA@@?=BABA@@AA@@@B@BA=A=AA@=AAAAAB<A?AA@B??>B?@@@AA?BBA?BA>?AB@AABA=AA?ABAA@AA>AA@AA@@B@A@ABAA	NM:i:5	MD:Z:0G0T0T0A3A142	MC:Z:150M	AS:i:144	XS:i:20	RG:Z:na12878_wgs"
   read = r"A00234:1477:HLFYYDSXC:2:2449:3152:5948	163	1	10001	2	9S110M16S	=	9998	119	CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACACTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCATAACCCGAACCG	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:,FF:,F:FF,F:,FF:FF,F,,:,,FF,,:,:,:,,	NM:i:1	MD:Z:28C81	MC:Z:12S111M1I11M	AS:i:105	XS:i:104	RG:Z:na12878-I2	XA:Z:1,+10028,81M1I35M18S,2;4,-191043982,16S67M1D44M3D8M,5;"

   t = time()
   r = __extract_reference_by_MD(read)
   print(r)
   print(read.split()[9])
   print(time() - t)


   read1 = r"00000000000:S:00000000000000:1:000001:R004:C065	99	1	15932376	60	150M	=	15932400	174	GGCACAAAGTCAATATTTGGAAATGTTAGCAATCATTTGGAAGGTGCTTGACGCACAGCGTTGATCGATAAATGGGACCTATTCTTCTTGGTAGCACTGCTGGTCTAGAGCAGGATTCTGGTTTCTGAGCTCTCCCTGNNCNNAACNTCC	AA@BAB@@A@AAAABAAAA@BAABA>AAAAAAAA@AABBA@A@A?@A?BAAABA@AAABAAA@AB@A@@@@AAAAA@AABAAAAAA@AABABA?AAAAAAABAABBAA?A?ABA@BAA@A@A@@AAB@A?AB@A=<@A!!2!!B;B!=:?	NM:i:5	MD:Z:138C0T1T0C3C3	MC:Z:150M	AS:i:140	XS:i:19	RG:Z:na12878_wgs"
   read2 = r"00000000000:S:00000000000000:1:000001:R004:C065	147	1	15932400	60	150M	=	15932376	-174	NNNNGCANTCATTTGGAAGGTGCTTGACGCACAGCGTTGATCGATAAATGGGACCTATTCTTCTTGGTAGCACTGCTGGTCTAGAGCAGGATTCTGGTTTCTGAGCTCTCCCTGCTCTCAACCTCCACTCCACCCCATCCTGACTTAAAA	!!!!:A@!@AA:A?>?@AB@?=@/@@BA>>AA=@BAA?A??@?BA7<?AAB??AB?AA@@?=BABA@@AA@@@B@BA=A=AA@=AAAAAB<A?AA@B??>B?@@@AA?BBA?BA>?AB@AABA=AA?ABAA@AA>AA@AA@@B@A@ABAA	NM:i:5	MD:Z:0G0T0T0A3A142	MC:Z:150M	AS:i:144	XS:i:20	RG:Z:na12878_wgs"
   t = time();r = get_base_info(read1, read2, 'test', 'test');print(len(r), time() - t)

   sam_file = '~/server/personal/projects/error_profile_all_in_one/evo_na12878/evo_202311_na12878_wgs_pe150_2023/sam/SalusEvo2_chrY_queryname.sam'
   for read1, read2 in iter_queryname_sam(sam_file):
      t = time()
      r = get_base_info(read1, read2, 'test', 'test', dummy_real_base = 'test')
      print(len(r), time() - t)
      i += 1
      if i % 100 == 0:
         break



   read1 = r"C2305160001:S:PRM32401030122:1:000000:R001:C019	83	12	133219225	60	146M4S	=	133219198	-173	GGAATTCCTCCAAGACAGGAATTTCACTGGCCAGCCTCTTCAGCTCCCAGCTGGACTGAACAGCGATGAGTGTGGGCCCCCGGCGCTCCTCCTGGGTCAAGGCAAAATGGAAGAAAANACCTGGGTGGACCCAGCCTCAAAGAAGACTTT	@A@AA@A2@B@A@@AA@AA?A@?<?BA@AAAAA?@A?A<@AA?AA@B0AAAAA@A@BA@@@AA@AAA=BA@A@AA=5@@@BA@@A@@AA@@A@@@@BBA@@A@BA@A?@B@A@@A@@!@A@@A@@@AAA@?B@@@AAB@ABB@AA@@;AA	NM:i:1	MD:Z:117G28	MC:Z:4S146M	AS:i:144	XS:i:0	RG:Z:sample2_salus_500M"
   read2 = r"C2305160001:S:PRM32401030122:1:000000:R001:C019	163	12	133219198	60	4S146M	=	133219225	173	ACTGTGTCAGCCACACAGATAGGCACCAGTGGGAATTCCTCCAAGACAGGAATTTCACTGGCCAGCCTCTTCAGCTCCCAGCTGGACTGAACAGCGATGAGTGTGGGCCCCCGGCGCTCCTCCTGGGTCAAGGCAAAATGGAAGAAAAGA	@@@@@A@@@AAA@@B@@@A@?A@B@A@AA@@@@BAA@A@A@AAAAB@@@A@@BAA@@BBAA@@@AB@A@A@@?@@B@@@A@B@AA?@AA@@BAB@BA@@??A@@AA@AA@@AA@AB@?@A@@BB@@@@@ABA@@?BBA@@?AAA:A@A@@	NM:i:0	MD:Z:146	MC:Z:146M4S	AS:i:146	XS:i:0	RG:Z:sample2_salus_500M"
   t = time();r = get_base_info(read1, read2, 'test', 'test');print(len(r), time() - t)
   # t = time();r = get_base_info(read2, read1, 'test', 'test');print(len(r), time() - t)

   i = 1
   '''
   '''
   chr_str = 'chr10'
   pos = 52976919
   bam_file = '/home/user/server/personal/projects/salus_analysis_20240619/KDTest22_60_meth_bwa/new/KDTest22.unmap.sort.filter.bam'
   bam_af = pysam.AlignmentFile(bam_file)
   base_real = __get_major_allele(bam_af, chr_str, pos, {}) # 可能为‘’，表示deletion

   i = 1
   '''

   read1 = r"00000000000:S:00000000000000:2:000542:R005:C056	99	ad36f57dd29d43c6_1	192962	60	150M	=	193068	256	ACGAAACCCAAATCGTCAACGTAATGACGATAAAAGCGCAGATCGCCAGAGTCAGCACGCGTGAGTTAATTTTCTTGGTATCAATAATTCGGCTTAACAGATTGAGAATAATGCCTTTAATGGCCTCGTGGAAACCGAGATAAATGCCAA	A@@AAAAAA@@AABA@AAB@AAAAA<@@B@AAA@A@A<BB<@@@AA@AB@?@AAAA@A?@BBB@>@B?ABAABAAAA2AA@A@@@A@@A@BAAAA@@@BBBAA@@@;AA@?A@ABA@?@ABB@@A??A<A@@A>B@@?B@AAA@2B@@AB	NM:i:0	MD:Z:150	MC:Z:150M	AS:i:150	XS:i:0	RG:Z:Ecoli_2"
   read2 = r"00000000000:S:00000000000000:2:000542:R005:C056	147	ad36f57dd29d43c6_1	193068	60	150M	=	192962	-256	AATCATGCCTTTAATGGCCTCGTGGAAACCGAGATAAATGCCAAAGAATGCGGTCAGTACGGCAAAGATATTAAGCACCGTAGAGGTGATATGAATGATATGCCCAGGGATCACCTGCGCGGCCAGCGCCAGTGCTGAGATATTTTGTTC	A4?4A@<@<@A@,B?@@A@AB@A@ABA@=@A@AAA@A@A@ABA@@ABA?@AA@A?@AAA?BA?@:B@BAA@AA@@@=7A>AA@@@@?ABA@@@<AA@@BAAAA@AB@@A?AA@?@@A@?A@A?A@@A@BB@B@ABBAA@A@ABA5@AA;A	NM:i:1	MD:Z:3A146	MC:Z:150M	AS:i:146	XS:i:0	RG:Z:Ecoli_2"
   t = time();r = get_base_info(read1, read2, 'test');print(len(r), time() - t)

   print(r)
