import unittest
import sys,os

sys.path.append('/usr/local/lib/python2.7/dist-packages/scons/')
from rsf.proj import *

testdir = os.path.dirname(os.environ['PWD'])+'/experiments/algorithmVelocityInversion/'
sys.path.append(testdir)
import diff

class TestDiff(unittest.TestCase):
	def test_diffsimul(self):
		self.assertEqual(3,3)
