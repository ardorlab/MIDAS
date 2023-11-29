# Unit test - MIDAS

import sys
sys.path.append('../')
sys.path.append('../../')

import unittest
from applications.ncsu_core import Extractor
test_file = open("test_file.out", "r")
file_lines = test_file.readlines()
test_file.close()

# Initiating variables to be tested.
keff_list = []
boron_list = []
EFPD_list = []
FDH_list = []
peak_list = []
exposure_list = []
boron_list = []

# Obtaining variables to be tested.
keff_list = Extractor.core_keff_list(file_lines)
EFPD_list = Extractor.efpd_list(file_lines)
boron_list = Extractor.boron_list(file_lines)
FDH_list = Extractor.FDH_list(file_lines)
peak_list = Extractor.pin_peaking_list(file_lines)
exposure_list = Extractor.burnup_list(file_lines)

# Known variables
known_keff_list = [1.0, 1.00001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.00001]
known_EFPD_list = [0.0, 13.0, 26.0, 52.0, 104.0, 156.0, 208.0, 260.0, 312.0, 364.0, 416.0, 459.4]
known_boron_list = [767.6, 742.2, 787.1, 854.8, 885.4, 838.1, 745.7, 621.4, 476.5, 317.8, 151.4, 10.0]
known_FDH_list = [1.414, 1.385, 1.369, 1.375, 1.333, 1.33, 1.325, 1.317, 1.305, 1.294, 1.283, 1.273]
known_peak_list = [1.783, 1.783, 1.785, 1.812, 1.677, 1.558, 1.532, 1.511, 1.466, 1.424, 1.422, 1.423]
known_exposure_list = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 17.669]

class TestExtractor(unittest.TestCase):
  
  def test_core_keff_list(self):
    self.assertEqual(keff_list, known_keff_list)
  
  def test_efpd_list(self):
    self.assertEqual(EFPD_list, known_EFPD_list)
    
  def test_boron_list(self):
    self.assertEqual(boron_list, known_boron_list)
    
  def test_FDH_list(self):
    self.assertEqual(FDH_list, known_FDH_list)
    
  def test_pin_peaking_list(self):
    self.assertEqual(peak_list, known_peak_list)
    
  def test_burnup_list(self):
    self.assertEqual(exposure_list, known_exposure_list)
    
if __name__ == '__main__':
  unittest.main()






