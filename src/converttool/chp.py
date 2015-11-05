#!/usr/bin/env python

import re
import sys


def replace(s,key):
	#s= re.sub(r'([ \+\*\/\-(),\=])'+key+r'([ \+\*\/\-(),$\=])',r"\1md%"+s+r"\2",s)
	#s= re.sub(r'^'+key+'([ \+\*/\-(),$=])',"md%uv\\1",s)
	#s= re.sub(r'([ +*/-(),=])'+key+r'([ +*/-(),=])',r"\1qmasmd%"+s+r"\2",s)
	#s= re.sub(r'([^\w)'+key+r'([^\w])',r"\1qmasmd%"+s+r"\2",s)
	s= re.sub('([^0-9A-Za-z_]|^)'+key+'([^0-9A-Za-z_]|$)',r'\1QMD%'+key+r'\2',s)
	s= re.sub('([^0-9A-Za-z_%]|^)'+key+'([^0-9A-Za-z_]|$)',r'\1QMD%'+key+r'\2',s)
	#s= re.sub('^'+key+'([^0-9A-Za-z_])',r'QMD%'+key+r'\1',s)
	return s

if __name__ == "__main__":

	filename=sys.argv[1]
	f=open(filename,'r')
	lines=f.readlines()
	f.close()

	wlist=['natom', 'nloopa', 'nloopc','loopa','loopc', 'imd', 'imdc', 'isd_cell', 'nkatm', 'ifrcf', 'iamax', 'uv', 'uvo', 'vuv', 'bv', 'omega', 'omegai', 'mconv', 'tote', 'toter0', 'strs', 'strso', 'extstrs', 'fth', 'fmax', 'sth', 'smax', 'tstep', 'mcell', 'guv', 'sdv_cell', 'duv0', 'lambda_cell', 'alphalm', 'iposfix', 'zatm', 'katm', 'ra', 'rr', 'rro', 'mass', 'mfac', 'frc', 'vrr', 'lgdiis', 'iuphess', 'igdiis', 'imod', 'imod0', 'imod2', 'imdgdiis', 'rah', 'gradh', 'hessi2', 'dampmu']

	for s in lines:	
		s=s.rstrip('\n')
		for w in wlist:
			s=replace(s,w)
		print s

