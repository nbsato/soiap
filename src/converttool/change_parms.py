#!/usr/bin/env python

def replace(lines):
	lines2=[]
	for s in lines:
		a=s.split()
		if len(a)==0:
			lines2.append(s)
			continue
		if a[0]=="parameter":
			lines2.append(s)
			lines2.append("type t_qmasmd")
			continue
		elif a[0]=="end" and a[1]=="module":
			lines2.append("end type")
			lines2.append("type(t_qmasmd):: QMD")
			lines2.append(s)
		else:
			lines2.append(s)	
	return lines2

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
        for s in lines:
		print s

