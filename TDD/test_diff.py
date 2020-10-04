import unittest
import sys,os

sys.path.append('/usr/local/lib/python2.7/dist-packages/scons/')
from rsf.proj import *

testdir = os.path.dirname(os.environ['PWD'])+'/experiments/algorithmVelocityInversion/'
sys.path.append(testdir)
from diff import diffsimul

class TestDiff(unittest.TestCase):
	def test_diffsimul(self):
		'''
		Function should raise a
		TypeError when velocities is
		not a tuple
		'''
		self.assertRaises(TypeError,
		diffsimul,'teste','testediff',1,1)

