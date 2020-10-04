# coding: utf-8
#
# inputScons.py (Madagascar Recipe)
# 
# Purpose: Functions to treat the input passed through
# command line in SCons.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo Dirack 04/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.


def castingParametersFromCmd(ARGLIST,param):
	'''
	casting of the parameters from
	command line
	:global ARGLIST: list of parameters passed through comand line in scons
	:param param: dictionary, parameters dictionary you want to cast
	'''

	for k, v in ARGLIST:
		if k in param.keys():
			if type(param[k]) is float:
				param[k] = float(v)
			elif type(param[k]) is int:
				param[k] = int(v)
			elif type(param[k]) is tuple:
				param[k] = eval(v)
			else:
				param[k] = v


