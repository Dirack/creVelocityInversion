# coding: utf-8
#
# diffSimulationAndMigration.py (Madagascar Recipe)
# 
# Purpose: Generate diffraction simulated hyperbolas
# in the stacked section using iterative picking. And
# migrate it using focusing measure. The input is a 
# stacked section with only diffraction hiperbolas.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo Dirack 01/03/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

from rsf.proj import *
import string, sys, os

frect=0

def diffmig(diffSimulatedSection,
            diffMigratedSection,
            v0,
            nv,
            dv,
            nx,
            padx,
            nt,
            tmin=0,
            tmax=10,
            rect1=10,
            rect2=10,
            srect1=1,
            srect2=3,
            vslope=None,
            units='Km',
            f1=1,
            j3=1,
            dx=1,
            x0=0,
            beg1=0,
            frect1=0,
            frect2=0,
            an=1,
            nout=2048,
            vx0=None):
	'''
	Diffraction migration by focusing measure
	:param diffSimulatedSection: RSF filename, stacked Section with simulated diffractions
	:output diffMigratedSection: RSF filename, basename of the outputs
	|-> diffMigratedSection+pik: time migrated velocity 
	|-> diffMigratedSection+dip: dominant slope section
	|-> diffMigratedSection+slc: slice of the velocity cube
	|-> diffMigratedSection+vlf: velocity cube focusing measure
	|-> diffMigratedSection+foc: focusing measure
	|-> diffMigratedSection+sem: semblance
	'''

	if diffSimulatedSection == diffMigratedSection:
		raise ValueError("Filenames should not be equal")
		return False

	if frect==0:
		frect1=2*rect1

	if frect==0:
		frect2=2*rect2/j3
	
	# dominant slope of stacked section
	dip = diffMigratedSection+'-dip'

	# Calculate dominante slope, rect parameter
	# are the smooth radius
	Flow(dip,diffSimulatedSection,'dip rect1=%d rect2=%d' % (rect1,rect2))

	# Velocity analisys by focusing
	# Stolt Migration. It needs the cosine
	# tranform of the data. sfvczo does
	# the velocity continuation process
	# migrate data multiple times with a
	# constant velocity
	velcon='''
	pad n2=%d beg1=%d | cosft sign2=1 | put o3=0 | 
	stolt vel=%g |
	vczo nv=%d dv=%g v0=%g |
	transp plane=23 |
	cosft sign2=-1 |
	window n2=%d f1=%d |
	put o2=%g |
	transp plane=23
	''' % (padx,beg1,v0,nv,dv,v0,nx,beg1,x0)

	# Velocity continuation of diffraction
	vlf=diffMigratedSection+'-vlf'
	Flow(vlf,diffSimulatedSection,velcon)

	Flow(vlf+'q',diffSimulatedSection,
	    '''
	    math output="input*input" | %s | clip2 lower=0
	    ''' % velcon)

	# j#n is the jump parameter of sfwindow
	if j3 > 1:
		focus = 'window j3=%d | ' % j3
	else:
		focus = ''

	# Focusing analisys
	# local varimax is max when the velocity
	# collapses a diffraction
	focus = focus + '''
	focus rect1=%d rect3=%d |
	math output="1/abs(input)" |
	cut max1=%g | cut min1=%g
	''' % (frect1,frect2,tmin,tmax)

	# focusing measure
	foc=diffMigratedSection+'-foc'
	Flow(foc,vlf,focus)

	# semblance is the local varimax measure
	# sfmul multplies data
	# sfdivn does smooth data division
	sem = diffMigratedSection+'-sem'
	Flow(sem,[vlf,vlf+'q'],
	    '''
	    mul $SOURCE |
	    divn den=${SOURCES[1]} rect1=%d rect3=%d
	    ''' % (rect1,rect2))
	
	# Data after picking
	pik=diffMigratedSection+'-pik' 

	# if vslope, use it to mute data
	if vslope:
		pick = 'mutter x0=%g v0=%g half=n | ' % (vx0,vslope)
	else:
		pick = ''

	pick = pick + 'scale axis=2 | pick an=%g rect1=%d rect2=%d | window' % (an,rect1,rect2)

	if j3 > 1:
		pick2 = pick + '''
		transp |
		spline n1=%d d1=%g o1=%g |
		transp
		''' % (nx,dx,x0)
	else:
		pick2 = pick

	Flow(pik,sem,pick2)

	# slice of picking cube
	slc = diffMigratedSection + '-slc' 
	Flow(slc,[vlf,pik],'slice pick=${SOURCES[1]}')
	return True

def diffsimul(
	section,
	diffSimulatedSection,	
	velocities,
	numberOfReflectors
	):
	'''
	Simulate difraction hyperbolas in the picked points in the stacked section
	:param section: RSF filename, stacked section
	:output diffSimulatedSection: RSF filename, diffraction simulated section
	:param velocities: tuple, velocity layers should equal to the number of reflectors
	:param numberOfReflectors: int, number of reflectors to iterative picking
	'''

	if type(velocities) is not tuple:
		raise TypeError("velocities parameter should be a tuple")
		return False

	if type(numberOfReflectors) is not int:
		raise TypeError("numberOfReflectors parameter should be int")
		return False

	if numberOfReflectors <= 0:
		raise ValueError("numberOfReflectors should be major than 0")
		return False

	if section == diffSimulatedSection:
		raise ValueError("filenames should not be equal")
		return False

	if len(velocities) < numberOfReflectors:
		raise ValueError("number of velocities should be equal or major than number of reflectors")
		return False
	
	reflectorsList = []
	# Iterative picking - Loop over reflectors
	for i in range(numberOfReflectors):

		reflectorPickedPoints = 'reflectorPickedPoints-%i.txt' % i
		t0sFile = 't0s-%i' % i
		t0sAscii = 't0s-%i.asc' %i
		m0sFile = 'm0s-%i' % i
		m0sAscii = 'm0s-%i.asc' % i
		returnedSection = 'returnedSection-%i' % i
		diffSection = 'diffSection-%i' % i

		# Reflector iterative Picking
		Flow(reflectorPickedPoints,section,
			'''
			ipick
			''')

		# Next step is done with 'ascFormat.sh' Shell Script
		# please check if the script have permissions to execute
		# Build t0 coordinates file (pass 1 to generate t0s file)
		Flow(t0sAscii,reflectorPickedPoints,
			'''
			./ascFormat.sh 1 %s
			''' % (t0sAscii))

		Flow(t0sFile,t0sAscii,'sfdd form=native')

		# Build m0 coordinates file (pass 2 to generate m0s file)
		Flow(m0sAscii,reflectorPickedPoints,
			'''
			./ascFormat.sh 2 %s
			''' % (m0sAscii))

		Flow(m0sFile,m0sAscii,'sfdd form=native')

		# Diffraction simulation in stacked section
		Flow([returnedSection,diffSection],[section,t0sFile,m0sFile],
			'''
			diffsim diff=${TARGETS[1]} aperture=1
			t0=${SOURCES[1]} m0=${SOURCES[2]} v=%g freq=10 verb=y
			''' % (velocities[i]))

		section = returnedSection
		reflectorsList.append(diffSection)

	Flow(diffSimulatedSection,reflectorsList,
		'''
		add ${SOURCES[1:%d]}
		scale=%s
		''' % (
		len(reflectorsList),
		','.join(['1' for i in range(len(reflectorsList))])))

	return True
