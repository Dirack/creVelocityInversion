#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# viewFoldGraph.py (Python)
# 
# Purpose: Plot NIP sources in depth model.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo A C Neves (Dirack) 27/09/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

import matplotlib.pyplot as plt

x=[]
y=[]
f =  open("nipsources.txt","r")
for line in f:
	ys,xs = line.split()
	x.append(float(xs))
	y.append(float(ys))

f.close()

xx=[]
yy=[]
f =  open("result.txt","r")
for line in f:
	ys,xs = line.split()
	xx.append(float(xs))
	yy.append(float(ys))

f.close()

z=[]
f = open("surfsources.txt","r")
for line in f:
	z.append(float(line))

f.close()

z0 = [0]*len(z)

m1=[]
m01=[]
i=0
f = open("modelin1.txt","r")
for line in f:
	m1.append(float(line))
	m01.append(0.01*i)
	i=i+1

f.close()

m2=[]
m02=[]
i=0
f = open("modelin2.txt","r")
for line in f:
	m2.append(float(line))
	m02.append(0.01*i)
	i=i+1

f.close()

plt.xlim(0, 6)
plt.ylim(3,0)
plt.plot(list(x),list(y),'*',label="NIP sources")
plt.plot(list(xx),list(yy),'.',label="NIP result")
plt.plot(list(z),list(z0),'v',ms=7,label="m0")
plt.plot(list(m01),list(m1),ls='-',label="First reflector")
plt.plot(list(m02),list(m2),ls='-',label="Second reflector")
plt.legend(loc="lower right")
plt.title("Inversion Result - NIP sources location",color="black",fontsize=16)
plt.xlabel("Position (Km)",fontsize=16)
plt.ylabel("Depth (Km)",fontsize=16)
plt.show()
