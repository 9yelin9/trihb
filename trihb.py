#!/usr/bin/env python3

import os
num_thread = 1
os.environ['OMP_NUM_THREADS'] = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import h5py
import time
import argparse
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-sg', '--showlog', nargs='+', help='<path_data> : ShowLog\n')
parser.add_argument('-sl', '--showlat', nargs='+', help='<path_data> : ShowLat\n')
parser.add_argument('-so', '--showobs', nargs='+', help='<dir_data_list(sep=:)> : ShowObs\n')
args = parser.parse_args()                                                                      

mk_list = ['s', 'o', '^']
ls_list = ['-', '--', ':']
plt.rcParams.update({'font.size': 20})

def RegExSubInt(val, string):
	return int(re.sub(val, '', re.search('%s\d+' % val, string).group()))
def RegExSubFloat(val, string):
	return float(re.sub(val, '', re.search('%s[-]?\d+[.]\d+' % val, string).group()))

def SaveFig(path, footer, fig):
	paths = path.split('/')
	dn, fn = '/'.join(paths[:-1]), paths[-1]
	dn, fn = re.sub('data', 'fig', dn), re.sub('[.]txt|[.]dat|[.]h5', '_%s.png' % footer, fn)
	os.makedirs(dn, exist_ok=True)
	fig_name = dn+'/'+fn
	fig.savefig(fig_name)
	print('Figure saved at %s' % fig_name)

def ShowLog(path_data):
	fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

	with open(path_data, 'r') as f: log = np.loadtxt(f, skiprows=1)
	itv = len(log)//1000
	itr = log[:, 0]
	cum_e = [e/i for i, e in zip(itr, np.cumsum(log[:, 1]))]

	e_min, e_max, e_fin = min(cum_e), max(cum_e), cum_e[-1] 
	v = (e_fin - e_min) / (e_max - e_min)
	v = v - 0.45 if e_fin > (e_min + e_max)/2 else v + 0.15
	cax = ax.inset_axes([0.62, v, 0.3, 0.3])

	ax.axhline(cum_e[-1], c='r', label=r'$\langle E \rangle=%f$' % (log[-1, 1]))
	cax.axhline(cum_e[-1], c='r')

	ax.plot(itr/itv, cum_e)
	cax.plot(itr[-100:]/itv, cum_e[-100:])

	cax.grid(True)
	cax.ticklabel_format(axis='x', useOffset=False, style='plain')
	cax.ticklabel_format(axis='y', useOffset=False, style='plain')

	ax.legend(fontsize='x-small')
	ax.grid(True)
	ax.set_xlabel(r'itr $\times10^{%d}$' % (np.log10(itv)))
	ax.set_ylabel(r'$\langle E \rangle$')
	path_data_s = path_data.split('/')
	ax.set_title('/'.join(path_data_s[:2])+'/\n'+'/'.join(path_data_s[2:]), fontsize='small')
	SaveFig(path_data, 'log', fig)
	plt.show()

def ShowLat(path_data):
	L = RegExSubInt('L', path_data)
	T = RegExSubFloat('T', path_data)

	with h5py.File(path_data, 'r') as f:
		data = f['lat'][()]
		mz, rho1, rho2, ozz = np.array(list(f['obs'][()][0]))
	site = data['site']
	angle = []
	for theta, phi in data['angle']:
		angle.append([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
	angle = np.array(angle)
	rho = (-2/(3*np.sqrt(3)*L*L)) * (rho1 + rho2/T)

	fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
	norm = plt.Normalize(-1, 1)

	qv = ax.quiver(site[:, 0], site[:, 1], angle[:, 0], angle[:, 1], angle[:, 2], norm=norm, edgecolors='k', cmap='bwr', pivot='mid', scale_units='inches', scale=1.5, lw=1.0, width=0.01, headlength=3, headaxislength=3)
	#sc = ax.scatter(site[:, 0], site[:, 1], c=angle[:, 2], norm=norm, cmap='bwr', edgecolors='k', s=200)

	cax = ax.inset_axes([0.92, 0.05, 0.02, 0.3])
	cb = plt.colorbar(qv, cax=cax, format='%d')
	cb.ax.set_title('$S_z$', pad=15)
	cb.set_ticks([-1, 1])

	ax.grid(True)
	ax.set_axisbelow(True)
	ax.set_aspect('equal')
	ax.set_xticks(np.unique(site[:, 0]), labels=[])
	ax.set_yticks(np.unique(site[:, 1]), labels=[])
	ax.set_xlim([-0.5, max(site[:, 0])+0.5])
	ax.set_ylim([-0.5, max(site[:, 1])+0.5])
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	ax.set_title(r'$m_z=%.4f$ $\rho_s=%.4f$ $O_{zz}=%.4f$' % (mz, rho, ozz), fontsize='small')
	SaveFig(path_data, 'lat', fig)
	plt.show()

def ShowObs(dir_data_list):
	dir_data_list = dir_data_list.split(':')
	fig, ax = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

	for i, dir_data in enumerate(dir_data_list):
		label = dir_data.split('/')[1]

		L, M = RegExSubInt('L', dir_data), RegExSubInt('M', dir_data)

		df = []
		for fn in ['%s/%s' % (dir_data, fn) for fn in os.listdir(dir_data) if re.search('T.+.h5', fn)]:
			h, D, T = RegExSubFloat('h', fn), RegExSubFloat('D', fn), RegExSubFloat('T', fn)
			with h5py.File(fn, 'r') as f:
				mz, rho1, rho2, ozz = np.array(list(f['obs'][()][0]))
			rho = (-2/(3*np.sqrt(3)*L*L)) * (rho1 + rho2/T)
			df.append([h, D, T, mz, rho, ozz])
		df = pd.DataFrame(df, columns=['h', 'D', 'T', 'mz', 'rho', 'ozz'])
		ax[0].plot(df['T'], df['rho'], marker=mk_list[i], ls=ls_list[i], fillstyle='none', label=label)
		ax[1].plot(df['T'], df['ozz'], marker=mk_list[i], ls=ls_list[i], fillstyle='none', label=label)
		print('\n'+label, df, sep='\n')
	print()

	for i, obs in enumerate(['rho', 'ozz']):
		ax[i].grid(True)
		ax[i].set_xlabel('T')
		ax[i].set_ylabel(obs)
	ax[0].set_ylim([-0.01, 0.31])
	ax[1].set_ylim([-0.1, 3.1])
	ax[1].legend(fontsize='x-small')
	fig_name = 'fig/' + '_'.join([dn.split('/')[1] for dn in dir_data_list]) + '_obs.png'
	fig.savefig(fig_name)
	print('Figure saved at %s' % fig_name)
	plt.show()

if args.showlog: ShowLog(*args.showlog)
elif args.showlat: ShowLat(*args.showlat)
elif args.showobs: ShowObs(*args.showobs)
else: parser.print_help()

