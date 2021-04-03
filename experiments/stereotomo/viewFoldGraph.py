#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# viewFoldGraph.py (Python)
# 
# Objetivo: Visualizar gráfico de cobertura dos CMP's.
# 
# Site: https://dirack.github.io
# 
# Versão 1.0
# 
# Programador: Rodolfo A C Neves (Dirack) 27/09/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# Licença: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

import matplotlib.pyplot as plt

f =  open("x.txt","r")
xstr = filter(None,f.read().split('\n'))
x = map(lambda i: float(i),xstr)

f =  open("y.txt","r")
ystr = filter(None,f.read().split('\n'))
y = map(lambda i: float(i),ystr)

f = open("z.txt","r")
zstr = filter(None,f.read().split('\n'))
z = map(lambda i: float(i),zstr)
z0 = [0]*26

#print(len(list(z)))
#print(len(list(z0)))

plt.xlim(0, 6)
plt.ylim(3,0)
plt.plot(list(x),list(y),'*')
plt.plot(list(z),list(z0),'v')
plt.show()
