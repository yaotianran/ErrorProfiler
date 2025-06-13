# 其他函数
import os
import os.path as path
from time import time
import pysam
import random


# 利用samtools过滤和排序sam或者bam文件
# 输入XXXX.bam/sam 文件，在当前目录生成XXXX_filtered.bam(.bai) 和 XXXX_queryname.sam文件
def preprocess(input_file: str, thread: int) -> int:

   '''
   利用samtools过滤和排序sam或者bam文件

   输入XXXX.bam/sam 文件，在当前目录生成XXXX_filtered.bam(.bai) 和 XXXX_queryname.sam文件
   '''

   temp_name = 'tmp_' + str(time()).split('.')[1] + str(random.randint(1, 9999999))  #  形如'tmp_484399657712'

   if not os.path.exists(temp_name):
      os.makedirs(temp_name)

   # input_file ~ /share/home/test.bam(sam)
   # basename ~ test
   basename = path.splitext(path.split(input_file)[1])[0]
   input_file_str = path.realpath(path.expanduser(input_file))

   try:
      # 输入文件转为clean bam
      print('cleaning bam ...')
      pysam.view('-@',
                 str(thread),
                 '-bh',
                 '-f',
                 '3',
                 '-F',
                 '3852',
                 '-o',
                 '{}/{}.bam'.format(temp_name, temp_name),
                 input_file_str,
                 catch_stdout = False
                 )

      # 将 clean bam排序
      print('sortinging bam by coordinates...')
      pysam.sort('-@',
                 str(thread),
                 '-o',
                 '{}/{}.sort.bam'.format(temp_name, temp_name),
                 '{}/{}.bam'.format(temp_name, temp_name),
                 catch_stdout = False
                 )

      pysam.index('-@',
                  str(thread),
                  '{}/{}.sort.bam'.format(temp_name, temp_name),
                  catch_stdout = False
                  )

      # 按queryname排序
      print('sortinging bam by querynames...')
      pysam.sort('-@',
                 str(thread),
                 '-n',
                 '-o',
                 '{}/{}.queryname.sam'.format(temp_name, temp_name),
                 '{}/{}.bam'.format(temp_name, temp_name),
                 catch_stdout = False)



   except pysam.SamtoolsError as ex:
      print(ex)
      return None, None

   try:
      os.remove('{}/{}.bam'.format(temp_name, temp_name))
   except Exception as ex:
      print(ex)
      return None, None

   return '{}/{}.sort.bam'.format(temp_name, temp_name), '{}/{}.queryname.sam'.format(temp_name, temp_name)


# index the bam
def index(bam_file: str, thread: int) -> int:
   pysam.index('-@',
               str(thread),
               bam_file,
               catch_stdout = False
               )

   return 0

def count_to_unit(count: int, digit = 2) -> str:
   '''
   将一个大于10_000整数转化为带有K，M，G，T，P，E的字符串
   例如 11230转化为 “11.23K”， 7486522转化为“7.486M”
   '''
   if count < 10000:
      return(str(count))


   unit_lst = ['K', 'M', 'G', 'T', 'P', 'E', 'KE', 'ME', 'GE', 'TE', 'PE', 'EE', '']
   is_converted = False
   for i in range(len(unit_lst)):
      if round(count / 1000 ** (i + 1), digit) < 1:
         unit = unit_lst[i - 1]
         num = count / 1000 ** i
         is_converted = True
         break

   if not is_converted:
      return str(count)

   # print(round(num, 2), unit)
   return str(round(num, digit)) + unit


if __name__ == '__main__':
   # r = preprocess('/home/user/sda1/报告模板位点信息产品手册/output.bam', thread = 4)
   # print(r)
   r = count_to_unit(99, digit=2)
   print(r)

   r = count_to_unit(9_999, digit=2)
   print(r)

   r = count_to_unit(10_999, digit=2)
   print(r)

   r = count_to_unit(54_644_168_735, digit=2)
   print(r)

   r = count_to_unit(365_254_644_168_735, digit=2)
   print(r)

   r = count_to_unit(10_365_254_644_168_735, digit=2)
   print(r)

   r = count_to_unit(1_102_000_365_254_644_168_735, digit=2)
   print(r)


