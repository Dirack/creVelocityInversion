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

z=[]
f = open("surfsources.txt","r")
for line in f:
	z.append(float(line))

f.close()

z0 = [0]*len(z)

m=[]
m0=[]
i=0
f = open("modelin.txt","r")
for line in f:
	m.append(float(line))
	m0.append(0.01*i)
	i=i+1

f.close()

plt.xlim(0, 6)
plt.ylim(3,0)
plt.plot(list(x),list(y),'*')
plt.plot(list(z),list(z0),'v')
plt.plot(list(m0),list(m),'-')
plt.show()
