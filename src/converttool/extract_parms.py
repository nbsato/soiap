#!/usr/bin/env python
import re
import sys

def replace(lines):
	start=0
	lines2=[]
	for s in lines:
		a=s.split()
		if len(a)==0:
			lines2.append(s)
			continue
		if a[0]=="module":
			continue
		if a[0]=="parameter":
			start=1
			continue
		if a[0]=="end" and a[1]=="module":
			start=0
			continue
		if a[0][0]=="!":
			continue
		if start==1:
			s2=re.sub("^.*::","",s)
			s3=""
			ix=0
			for i in range(len(s2)):
				if s2[i]=="(":
					ix=1
				elif s2[i]==")":
					ix=0
					continue
				else:
					if ix==0:
						s3+=s2[i]
			lines2.append(s3)
	w=[]
	for s in lines2:
		s2=s.split(",")
		for s3 in s2:
			s4=s3.split()
			for s5 in s4:
				w.append(s5)
	print w

if __name__ == "__main__":

	filename="parms.f90"
        f=open(filename,'r')
        lines=f.readlines()
        f.close()
	lines2=[]
        for s in lines: 
                s=s.rstrip('\n')
		lines2.append(s)
	lines=replace(lines2)
	sys.exit(10)
        for s in lines:
		print s

