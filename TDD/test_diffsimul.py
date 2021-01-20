import unittest
import sys,os

sys.path.append('/usr/lib/scons/')
from rsf.proj import *

testdir = os.path.dirname(os.environ['PWD'])+'/experiments/algorithmVelocityInversion/'
sys.path.append(testdir)
from diffSimulationAndMigration import diffsimul

class TestDiffSimul(unittest.TestCase):

	def setUp(self):
		self.section = 'stackedSection'
		self.diffSection = 'diffSimulatedsection' 
		self.velocities = (1.5,1.6,2.0)
		self.numberOfReflectors = 3

	def test_diffsimul_velocities(self):
		'''
		Function should raise a
		TypeError when velocities is
		not a tuple
		'''
		self.assertRaises(TypeError,
		diffsimul,self.section,self.diffSection,1,self.numberOfReflectors)
	
	def test_diffsimul_numberOfReflectors_int(self):
		'''
		Function should raise a
		TypeError when numberOfReflectors
		is not an int
		'''
		self.assertRaises(TypeError,
		diffsimul,self.section,self.diffSection,self.velocities,2.0)

	def test_diffsimul_numberOfReflectors_positive(self):
		'''
		Function should raise a ValueError
		if numberOfReflectors is negative
		'''
		self.assertRaises(ValueError,
		diffsimul,self.section,self.diffSection,self.velocities,-1)

	def test_diffsimul_numberOfReflectors_major_0(self):
		'''
		Function should raise a ValueError
		if numberOfReflectors is not
		major than 0
		'''
		self.assertRaises(ValueError,
		diffsimul,self.section,self.diffSection,self.velocities,0)

	def test_diffsimul_filenames_not_equal(self):
		'''
		input filenames should not be equal
		'''
		self.assertRaises(ValueError,
		diffsimul,self.section,self.section,self.velocities,self.numberOfReflectors)

	def test_diffsimul_ok(self):
		'''
		Function normally returns True
		'''
		self.assertTrue(diffsimul(
		self.section,self.diffSection,self.velocities,1))
	
	def test_diffsimul_len_velocities(self):
		'''
		Number of velocities should be
		equal or major than 
		numberOfReflectors
		'''
		self.assertRaises(ValueError,
		diffsimul,self.section,self.diffSection,self.velocities,4)
