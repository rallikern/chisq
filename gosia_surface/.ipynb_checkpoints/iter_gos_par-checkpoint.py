"""
Nice Script for performing Chi square distribution with Gosia
"""

import sys
sys.path.insert(0, "/home/ralli/Physik/hieisolde/off_ana/mathematica_ana/")
import os
import shutil
import glob
from multiprocessing import Process
import math as mt
import pandas as pd
import numpy as np
import trans_strength as ts
import tran_prob as tp

OUT_NAME = 'UPL_me2_4_174'
DIR1 = '/home/ralli/Physik/hieisolde/off_ana/mathematica_ana/chisq/gosia_surface/input_files/'
DIR2A = DIR1
DIR2B = DIR1
DIR3 = DIR1
DIR4 = '/home/ralli/Physik/hieisolde/gosia_isolde/'
DIR5 = '/home/ralli/Physik/hieisolde/off_ana/mathematica_ana/chisq/gosia_surface/output_files/'
DIR6 = '/home/ralli/Physik/hieisolde/off_ana/mathematica_ana/chisq/gosia_surface/gosia_files/'
HEADER_TEMPL = DIR1 + 'gosia_header_dens_3'
HEADER = 'gosia_header_1'
DET_INP = 'miniball_array.inp'
YLD = '140nd_208pb.yld'
NTAP3 = pd.read_csv(DIR1+'ntap3.dat', delimiter="\t", skip_blank_lines=False, header=None)
NTAP4 = pd.read_csv(DIR1+'ntap4.dat', delimiter="\t", skip_blank_lines=False, header=None)
PART1 = pd.read_csv(DIR2A+'gosia_part1', delimiter="\t", skip_blank_lines=False, header=None)
PART2 = pd.read_csv(DIR2B+'gosia_part2_tot', delimiter="\t", skip_blank_lines=False, header=None)
CORR = pd.read_csv(DIR1+'gosia_corr', delimiter="\t", skip_blank_lines=False, header=None)
MAPI = pd.read_csv(DIR1+'gosia_map', delimiter="\t", skip_blank_lines=False, header=None)
REST = pd.read_csv(DIR1+'gosia_rest', delimiter="\t", skip_blank_lines=False, header=None)
MINI = pd.read_csv(DIR3+'gosia_mini', delimiter="\t", skip_blank_lines=False, header=None)
#TRAN_PAR = [140,-0.19,2.332,1.558,100,54]
TRAN_PAR = [140,-0.19,2.332,1.558,100,54]

LINES = []
LINES.append(tp.file_checker_min('parameter0,parameter0,parameter0', HEADER_TEMPL))
LINES.append(tp.file_checker_min('parameter1,parameter1,parameter1', HEADER_TEMPL))
LINES.append(tp.file_checker_min('parameter2,parameter2,parameter2', HEADER_TEMPL))
LINES.append(tp.file_checker_min('parameter3,parameter3,parameter3', HEADER_TEMPL))

LINES = np.array(LINES) - np.ones(len(LINES))


def minimizer(me, mx, direct, iterations=3):
    """
    Function for minimization
    """
    shutil.copyfile(HEADER_TEMPL, './'+HEADER)
    me2 = np.sign(TRAN_PAR[1])*ts.arachnomxe2(*TRAN_PAR[:2], 0, mx, 0, *TRAN_PAR[2:5], 0, TRAN_PAR[5], 0)[0]
    mm1 = mt.sqrt(5*ts.arachnobm1(*TRAN_PAR[:2], 0, mx, 0, *TRAN_PAR[2:5], 0, TRAN_PAR[5], 0)[0])
    #print(me2, mm1)
    values = []
    values.append('1,10'+3*(','+str(mx)))
    values.append('10,10'+3*(','+str(me)))
    values.append('2,10'+3*(','+str(me2)))
    values.append('2,10'+3*(','+str(mm1)))

    sts = os.system(DIR4+"gosia" + " < " + DET_INP)
    if sts == 2:
        os.chdir('../')
        shutil.rmtree('./'+direct, ignore_errors=True)
        return

    inp_ra = pd.read_csv(HEADER, delimiter="\t", skip_blank_lines=False, dtype=str, header=None)
    for lines in enumerate(LINES):
        inp_ra.loc[lines[1]] = values[lines[0]]
    inp_fi = pd.concat([inp_ra, NTAP3, PART1, PART2, CORR], axis=0, ignore_index=True)
    
    inp_fi.to_csv('gosia_nd.inp', index=False, sep=';', header=False, na_rep=' ')
    
    sts = os.system(DIR4+"gosia" + " < gosia_nd.inp")
    if sts == 2:
        os.chdir('../')
        shutil.rmtree('./'+direct, ignore_errors=True)
        return

    inp_fi = pd.concat([inp_ra, NTAP3, PART1, PART2, MAPI], axis=0, ignore_index=True)
    inp_fi.to_csv('gosia_nd.inp', index=False, sep=';', header=False, na_rep=' ')

    sts = os.system(DIR4+"gosia" + " < gosia_nd.inp")
    if sts == 2:
        os.chdir('../')
        shutil.rmtree('./'+direct, ignore_errors=True)
        return

    inp_fi = pd.concat([inp_ra, NTAP4, PART1, PART2, MINI], axis=0, ignore_index=True)
    inp_fi.to_csv('gosia_nd.inp', index=False, sep=';', header=False, na_rep=' ')

    sts = os.system(DIR4+"gosia" + " < gosia_nd.inp")

    if sts == 2:
        os.chdir('../')
        shutil.rmtree('./'+direct, ignore_errors=True)
        return

    print('successfull minimization')

    for _ in range(iterations):
        inp_fi = pd.concat([inp_ra, NTAP3, PART1, REST, PART2, CORR], axis=0, ignore_index=True)
        inp_fi.to_csv('gosia_nd.inp', index=False, sep=';', header=False, na_rep=' ')
        sts = os.system(DIR4+"gosia" + " < gosia_nd.inp")
        if sts == 2:
            os.chdir('../')
            shutil.rmtree('./'+direct, ignore_errors=True)
            break

        inp_fi = pd.concat([inp_ra, NTAP3, PART1, REST, PART2, MAPI], axis=0, ignore_index=True)
        inp_fi.to_csv('gosia_nd.inp', index=False, sep=';', header=False, na_rep=' ')
        sts = os.system(DIR4+"gosia" + " < gosia_nd.inp")
        if sts == 2:
            os.chdir('../')
            shutil.rmtree('./'+direct, ignore_errors=True)
            break

        inp_fi = pd.concat([inp_ra, NTAP4, PART1, REST, PART2, MINI], axis=0, ignore_index=True)
        inp_fi.to_csv('gosia_nd.inp', index=False, sep=';', header=False, na_rep=' ')
        sts = os.system(DIR4+"gosia" + " < gosia_nd.inp")
        if sts == 2:
            os.chdir('../')
            shutil.rmtree('./'+direct, ignore_errors=True)
            break

        print('successfull minimization')

    if sts == 2:
        os.chdir('../')
        shutil.rmtree('./'+direct, ignore_errors=True)
        return

    value1, value2 = str("{:0.3f}".format(me)), str("{:0.3f}".format(mx))
    if value1[0] == '0':
        value1 = value1[1:]
    if value2[0] == '0':
        value2 = value2[1:]

    val = value1+'.'+value2
    shutil.copyfile('140nd_208pb.out', DIR5+OUT_NAME+'/out_'+val)
    os.chdir('../')
    shutil.rmtree('./'+direct, ignore_errors=True)
    return

if __name__ == "__main__":  # confirms that the code is under main function
    try:
        os.mkdir(DIR5+OUT_NAME)
    except:
        print('Directory already exists!')

    for me1 in np.linspace(-0.60, -0.63, 4):
        mxs = np.linspace(0.10, 0.13, 4)
        me2_temp = me1
        directory = 'bashi_'+str(me1)
        def minimizer1(mx1):
            """
            Function for minimization
            """
            directory1 = directory+'_'+str(mx1)+'/'
            try:
                os.mkdir(directory1)
            except:
                shutil.rmtree(directory1, ignore_errors=True)
                os.mkdir(directory1)
            #fort = glob.glob(DIR6+'fort*')
            #for i in fort:
            #    i_final = i[len(DIR6):]
            #    shutil.copyfile(i, directory1 + i_final)
            #output = glob.glob(DIR6+'140nd_208pb.*')
            #for i in output:
            #   i_final = i[len(DIR6):]
            #shutil.copyfile(i, directory1 + i_final)
            directory_gos = 'gosia_files/'
            for i in [YLD, DET_INP]:
                shutil.copyfile(directory_gos + i, directory1 + i)
            os.chdir('./'+directory1)
            return minimizer(me2_temp, mx1, directory1, iterations=3)

        procs = []
        for mx2 in mxs:
            # print(name)
            proc = Process(target=minimizer1, args=(mx2,))
            procs.append(proc)
            proc.start()

        # complete the processes
        for proc in procs:
            proc.join()
