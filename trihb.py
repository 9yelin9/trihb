#!/usr/bin/env python3

import os
num_thread = 1
os.environ['OMP_NUM_THREADS'] = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import time
import argparse
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-gs', '--genscript', nargs='+', type=float, help='<h> <D> <T> <v_init> <v_final> <dv> : GenScript\n* Enter -1 for the variable to set the range\n')
parser.add_argument('-gj', '--genjob', nargs='+', help='<dir_data> <init_mode> <Meq> <Mmc> <L> [q=openmp.q@phase09]: GenJob\n')
parser.add_argument('-so', '--showobs', nargs='+', help='<dir_data_list(sep=,)> <x=h,D,T> <y=mz,rho,ozz> : ShowObs\n')
args = parser.parse_args()                                                                     

def GenScript(h, D, T, v_init, v_final, dv):
	if v_init > v_final: dv = -dv
	v_range = np.arange(v_init, v_final+dv, dv)
	v_target = ''
	dir_data = ''
	
	v_arr = []
	for v, v_str in zip([h, D, T], ['h', 'D', 'T']):
		if v < 0:
			v_arr.append(v_range)
			v_target = v_str
		else:
			v_arr.append([v for _ in v_range])
			dir_data += '%s%.4f_' % (v_str, v)
	v_arr = np.array(v_arr).T
	
	dir_data += '%sF%.4f_%sT%.4f' % (v_target, v_init, v_target, v_final)
	tail = len([d for d in os.listdir('data/') if re.match(dir_data, d)])
	if tail: dir_data += '_%d' % tail
	os.makedirs('data/'+dir_data, exist_ok=True)
		
	fn = 'data/%s/script.txt' % dir_data
	with open(fn, 'w') as f:
		f.write('%d\n' % len(v_arr))
		np.savetxt(f, v_arr, fmt='%10f')

	print('GenScript(%s)' % fn)

def RegExSub(val, string):
	return float(re.sub(val, '', re.search('%s[-]?\d+[.]\d+' % val, string).group()))

def GenJob(dir_data, init_mode, Meq, Mmc, L, q='openmp.q@phase09'):
	Meq, Mmc, L = float(Meq), float(Mmc), int(L)

	dir_save = '%s/L%d_Mmc%.e/' % (dir_data, L, Mmc) if Meq == 0 else '%s/L%d_Meq%.e_Mmc%.e/' % (dir_data, L, Meq, Mmc)
	os.makedirs(dir_save, exist_ok=True)
	os.makedirs(dir_save + 'job/', exist_ok=True)
	os.makedirs(dir_save + 'log/', exist_ok=True)

	with open('%s/script.txt' % dir_data, 'r') as fs, open('job.sh', 'r') as fj:
		for line_s in fs:
			fn = '%s/job/%s.sh' % (dir_save, '_'.join(['%s%.4f' % (a, float(v)) for a, v in zip(['h', 'D', 'T'], line_s.split())]))
			with open(fn, 'w') as f:
				for line_j in fj:
					if re.search('#[$] -q', line_j): f.write('#$ -q %s\n' % q)
					elif re.search('#[$] -o', line_j): f.write('#$ -o %s/log/$JOB_NAME.log\n' % dir_save)
					elif re.search('###', line_j): f.write('./mc %s %s %.e %.e %d %s' % (dir_data, init_mode, Meq, Mmc, L, line_s))
					else: f.write(line_j)
				fj.seek(0)

	print('GenJob(%s)' % (dir_save))

def ShowObs(dir_data_list, x, y):
	dir_data_list = dir_data_list.split(',')

	mk_list = ['s', 'o', '^']
	ls_list = ['-', '--', ':']
	fig, ax = plt.subplots(figsize=(8, 6))

	for dn, mk, ls in zip(dir_data_list, mk_list, ls_list):
		label = dn.split('/')[-2]
		fn = dn + '/obs.txt'
		if not os.path.isfile(fn): os.system('./obs %s' % dn)
		df = pd.read_csv(fn, sep='\s+', comment='#', names=['h', 'D', 'T', 'mz', 'rho', 'ozz'])
		ax.plot(df[x], df[y], marker=mk, ls=ls, fillstyle='none', label=label)

	ax.grid(True)
	ax.legend()
	ax.set_title(re.sub('_vs.+', '', dir_data_list[0].split('/')[1]))
	ax.set_xlabel(x)
	ax.set_ylabel(y)
	plt.show()

if args.genscript: GenScript(*args.genscript)
elif args.genjob: GenJob(*args.genjob)
elif args.showobs: ShowObs(*args.showobs)
else: parser.print_help()
