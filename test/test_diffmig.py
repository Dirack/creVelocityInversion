import unittest
import sys,os

sys.path.append('/usr/lib/scons/')
from rsf.proj import *

testdir = os.path.dirname(os.environ['PWD'])+'/experiments/algorithmVelocityInversion/'
sys.path.append(testdir)
from diffSimulationAndMigration import diffmig

class TestDiffMig(unittest.TestCase):

	def setUp(self):
		self.diffSimulatedSection = 'diffSimulatedSection'
		self.diffMigratedSection = 'diffMigratedsection' 
		self.param = {'diffSimulatedSection': self.diffSimulatedSection,
		'diffMigratedSection': self.diffMigratedSection,
		'v0': 1.4,
		'nv':100,
		'dv':0.01,
		'nx':201,
		'padx':1000,
		'nt':1001,
		'tmin':0,
		'tmax':4,
		'rect1':10,
		'rect2':10,
		'srect1':1,
		'srect2':3,
		'vslope':None,
		'units':'Km',
		'f1':1,
		'j3':1,
		'dx':0.025,
		'x0':0,
		'beg1':0,
		'frect1':0,
		'frect2':0,
		'an':1,
		'nout':2048,
		'vx0':None}

		self.tmp={}

	def test_diffmig_ok(self):
		'''
		Function normally returns True
		if succed
		'''
		self.assertTrue(
		diffmig(**(self.param)))

	def test_diffmig_filenames(self):
		'''
		Filenames should not be equal
		function should raise a 
		ValueError
		'''
		self.tmp = self.param
		self.tmp['diffMigratedSection'] = self.param['diffSimulatedSection']
		self.assertRaises(ValueError,
		diffmig,**self.tmp)
