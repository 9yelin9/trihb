#!/usr/bin/env python3

import os
num_thread = 1
os.environ['OMP_NUM_THREADS'] = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import time
import numpy as np
import pandas as pd 
import inspect
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-s', '--script', nargs='+', type=float, help='<h> <D> <T> <v_start> <v_end> <dv> : GenScript\n* Enter -1 for the variable to set the range\n')
parser.add_argument('-j', '--job', nargs='+', help='<dir_init> <init_mode> <dir_data> <L> <M> [q=openmp.q@phase09]: GenJob\n')
parser.add_argument('-o', '--obs', nargs='+', help='<dir_data_list(sep=,)> <x=h,D,T> <y=mz,rho,ozz> : ShowObs\n')
args = parser.parse_args()                                                                     

def GenScript(h, D, T, v_start, v_end, dv):
	v_range = np.arange(v_start, v_end+dv, dv)
	v_list = []
	for v in [h, D, T]:
		if v < 0: v_list.append(v_range)
		else: v_list.append([v for _ in v_range])
	v_list = np.array(v_list).T
	
	dir_data = 'data/h%.4f_D%.4f_T%.4f_vs%.4f_ve%.4f_dv%.4f' % (h, D, T, v_start, v_end, dv)
	fn = '%s/script.txt' % dir_data

	os.makedirs(dir_data, exist_ok=True)
	np.savetxt(fn, v_list, fmt='%10f')

	print('GenScript(%s)' % fn)

def GetVal(val, string):
	return float(re.sub(val, '', re.search('%s[-]?\d+[.]\d+' % val, string).group()))

def GenJob(dir_init, init_mode, dir_data, L, M, q='openmp.q@phase09'):
	arg, _, _, val = inspect.getargvalues(inspect.currentframe())
	tm = time.strftime('%m%d%H%M', time.localtime(time.time()))

	dir_save = '%s/L%s_M%s_tm%s' % (dir_data, L, M, tm)
	os.makedirs(dir_save, exist_ok=True)

	info = '%s/info.txt' % dir_save
	with open(info, 'w') as f:
		for a in arg[:-1]: f.write('%s %s\n' % (a, val[a]))
		f.write('%s %s\n' % ('tm', tm))

	"""
	cnt = 0
	for v in ['h', 'D', 'T']:
		if GetVal(v, dir_data) < 0:
			for fn in os.listdir('job/'):
				if re.search(dn+v, 'job/'+fn): cnt += 1
			dn += '%s%d_' % (v, cnt)
			break
	dn += 'L%s_M%s' % (L, M)
	"""
	dn = 'job/%s' % (re.sub('data_', '', re.sub('__', '_', re.sub('/', '_', dir_save))))
	os.makedirs(dn, exist_ok=True)	
	os.makedirs(re.sub('job', 'log', dn), exist_ok=True)

	with open('%s/script.txt' % dir_data, 'r') as fs, open('job/default.sh', 'r') as fd:
		for line0 in fs:
			fn = '%s/%s.sh' % (dn, '_'.join(['%s%.4f' % (a, float(v)) for a, v in zip(['h', 'D', 'T'], line0.split())]))
			with open(fn, 'w') as f:
				for line in fd:
					if re.search('#[$] -q', line): f.write('#$ -q %s\n' % q)
					elif re.search('#[$] -o', line): f.write('#$ -o %s/$JOB_NAME.log\n' % re.sub('job', 'log', dn))
					elif re.search('###', line): f.write('./mc %s %s %s %s %s %s' % (dir_init, init_mode, dir_save, L, M, line0))
					else: f.write(line)
			fd.seek(0)

	print('GenJob(%s, %s)' % (dir_save, dn))

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

if args.script: GenScript(*args.script)
elif args.job: GenJob(*args.job)
elif args.obs: ShowObs(*args.obs)
else: parser.print_help()
