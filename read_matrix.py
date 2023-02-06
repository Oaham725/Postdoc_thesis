"""
Author: Hao Ma
Last update: 2022-09-23
Mail: oaham@xmu.edu.cn
Usage: 计算预共振拉曼的极化率张量导数
"""
import re
from math import sqrt
import numpy as np
import pandas as pd
import time
import glob
import os
import linecache

def noE(line):
    #将Fortran输入的值改成python可以读写的形式，正则表达式是最好用的。
    line = re.sub(r'D', "E", line)
    Newline = re.findall('([-]?\d+[.]\d+[E][+-]\d+)',line)
    return Newline

'''需要处理的文件路径，目前需要将对应的
Property number 1 -- Alpha(-w,w) derivatives, frequency  2    0.071980:
后续文件另存为一个新的文件。读写顺序为横向拼接，每3x3为一个极化率张量矩阵，共有63个。
'''
Files = glob.glob(r"C:\Users\Administrator\Desktop\matrix.txt")

for file in Files:
    f = open(file)
    path = os.path.realpath(file)
    parent = os.path.dirname(os.path.realpath(file))
    portion = os.path.splitext(file)
    newname = portion[0] + ".pol"
    lines = f.readlines()
    i = int(len(lines)/4)
    line1 = []
    line2 = []
    line3 = []
    for j in range(0, i):
        line1 = line1 + noE(lines[4 * j + 1])
        # print(line1)
        line2 = line2 + noE(lines[4 * j + 2])
        line3 = line3 + noE(lines[4 * j + 3])
    res1 = list(map(float, line1))
    # print(res1)
    res2 = list(map(float, line2))
    res3 = list(map(float, line3))
    with open(newname, "w") as aha:
        aha.write("POLA" + "\n")
        for num in range(63): # 更改原子数乘以三，对应横向维度
            alpha1 = '%.5f' % float(res1[3*num])
            alpha2 = '%.5f' % float(res2[3*num])
            alpha3 = '%.5f' % float(res2[3 * num+1])
            alpha4 = '%.5f' % float(res3[3 * num])
            alpha5 = '%.5f' % float(res3[3 * num+1])
            alpha6 = '%.5f' % float(res3[3 * num+2])
            #对整体数值进行调整，因为fortran的读写要求是：数值（6F10.6）格式为（20A4）
            ''' xx '''
            if abs(float(alpha1)) > 100:
                alpha1 = " " + '%.4f' % float(res1[3*num])
            if float(alpha1) > 0:
                alpha1 = " " + alpha1
            if abs(float(alpha1)) > 10 and abs(float(alpha1)) < 100:
                alpha1 = " " + alpha1
            if abs(float(alpha1)) < 10:
                alpha1 = "  " + alpha1
            ''' xy '''
            if abs(float(alpha2)) > 100:
                alpha2 = " " + '%.4f' % float(res2[3 * num])
            if float(alpha2) > 0:
                alpha2 = " " + alpha2
            if abs(float(alpha2)) > 10 and abs(float(alpha2)) < 100:
                alpha2 = " " + alpha2
            if abs(float(alpha2)) < 10:
                alpha2 = "  " + alpha2
            ''' yy '''
            if abs(float(alpha3)) > 100:
                alpha3 = " " + '%.4f' % float(res2[3 * num + 1])
            if float(alpha3) > 0:
                alpha3 = " " + alpha3
            if abs(float(alpha3)) > 10 and abs(float(alpha3)) < 100:
                alpha3 = " " + alpha3
            if abs(float(alpha3)) < 10:
                alpha3 = "  " + alpha3
            ''' xz '''
            if abs(float(alpha4)) > 100:
                alpha4 = " " + '%.4f' % float(res3[3 * num])
            if float(alpha4) > 0:
                alpha4 = " " + alpha4
            if abs(float(alpha4)) > 10 and abs(float(alpha4)) < 100:
                alpha4 = " " + alpha4
            if abs(float(alpha4)) < 10:
                alpha4 = "  " + alpha4
            ''' yz '''
            if abs(float(alpha5)) > 100:
                alpha5 = " " + '%.4f' % float(res3[3 * num + 1])
            if float(alpha5) > 0:
                alpha5 = " " + alpha5
            if abs(float(alpha5)) > 10 and abs(float(alpha5)) < 100:
                alpha5 = " " + alpha5
            if abs(float(alpha5)) < 10:
                alpha5 = "  " + alpha5
            ''' zz '''
            if abs(float(alpha6)) > 100:
                alpha6 = " " + '%.4f' % float(res3[3 * num + 2])
            if float(alpha6) > 0:
                alpha6 = " " + alpha6
            if abs(float(alpha6)) > 10 and abs(float(alpha6)) < 100:
                alpha6 = " " + alpha6
            if abs(float(alpha6)) < 10:
                alpha6 = "  " + alpha6
            ''' output '''
            aha.write( alpha1  + alpha2  + alpha3  + alpha4  + alpha5  + alpha6 + "\n")
        aha.write("\n")
