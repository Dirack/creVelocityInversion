import unittest
import sys,os

sys.path.append('/usr/local/lib/python2.7/dist-packages/scons/')
from rsf.proj import *

testdir = os.path.dirname(os.environ['PWD'])+'/experiments/algorithmVelocityInversion/'
sys.path.append(testdir)
from diffSimulationAndMigration import diffmig

class TestDiffMig(unittest.TestCase):

	def setUp(self):
		self.diffSimulatedSection = 'diffSimulatedSection'
		self.diffMigratedSection = 'diffMigratedsection' 

	def test_diffmig_ok(self):
		'''
		Function normally returns True
		'''
		self.assertTrue(
		diffmig(self.diffSimulatedSection,
		self.diffMigratedSection,
		v0=1.4,
		nv=100,
		dv=0.01,
		nx=201,
		padx=1000,
		nt=1001,
		tmin=0,
		tmax=4,
		rect1=10,
		rect2=10,
		srect1=1,
		srect2=3,
		vslope=None,
		units='Km',
		f1=1,
		j3=1,
		dx=0.025,
		x0=0,
		beg1=0,
		frect1=0,
		frect2=0,
		an=1,
		nout=2048,
		vx0=None))

	
