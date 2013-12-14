#!/usr/bin/env python
#-*- coding:utf-8 -*-

import unittest
import sys
sys.path.append( "../" )
from demultadapt import *

class TestLevenshtein_selector(unittest.TestCase):

    def test_single(self):
        lsof = Levenshtein_selector( [ ("ATCGCA", 0),
                                                 ("CCAGTG", 1),
                                                 ("GGTAAT", 2), ], True, 0.75)
        
        self.assertEqual( lsof.select( "CCAGTG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "ATCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GGTAAT" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "CCAGGG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "TTCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GCTAAT" ), ("GGTAAT", 2) )

        self.assertIs( lsof.select( "AAAAAA" ), None )
        self.assertIs( lsof.select( "CCAGCA" ), None )

    def test_paired(self):
        lsof = Levenshtein_selector( [ ("ATCGCA", 0),
                                       ("CCAGTG", 1),
                                       ("GGTAAT", 2), ], False, 0.75)
        
        self.assertEqual( lsof.select( "CCAGTG", "CCAGTG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "ATCGCA", "ATCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GGTAAT", "GGTAAT" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "CCAGGG", "CCAGGG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "TTCGCA", "TTCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GCTAAT", "GCTAAT" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "CCAGGG", "AAAAAA" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "TTCGCA", "AAAAAA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GCTAAT", "AAAAAA" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "AAAAAA", "CCAGTG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "AAAAAA", "ATCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "AAAAAA", "GGTAAT" ), ("GGTAAT", 2) )

        self.assertIs( lsof.select( "AAAAAA", "AAAAAA" ), None )
        self.assertIs( lsof.select( "CCAGTG", "ATCGCA" ), None )
        
        
class TestStd_selector(unittest.TestCase):

    def test_single(self):
        lsof = Std_selector( [ ("ATCGCA", 0),
                               ("CCAGTG", 1),
                               ("GGTAAT", 2), ], True )
        
        self.assertEqual( lsof.select( "CCAGTG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "ATCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GGTAAT" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "CCAGGG" ), None )
        self.assertEqual( lsof.select( "TTCGCA" ), None )
        self.assertEqual( lsof.select( "GCTAAT" ), None )


    def test_paired(self):
        lsof = Std_selector( [ ("ATCGCA", 0),
                               ("CCAGTG", 1),
                               ("GGTAAT", 2), ], False )
        
        self.assertEqual( lsof.select( "CCAGTG", "CCAGTG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "ATCGCA", "ATCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GGTAAT", "GGTAAT" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "CCAGGG", "CCAGGG" ), None )
        self.assertEqual( lsof.select( "TTCGCA", "TTCGCA" ), None )
        self.assertEqual( lsof.select( "GCTAAT", "GCTAAT" ), None )

        self.assertEqual( lsof.select( "CCAGTG", "AAAAAA" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "ATCGCA", "AAAAAA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "GGTAAT", "AAAAAA" ), ("GGTAAT", 2) )
        
        self.assertEqual( lsof.select( "AAAAAA", "CCAGTG" ), ("CCAGTG", 1) )
        self.assertEqual( lsof.select( "AAAAAA", "ATCGCA" ), ("ATCGCA", 0) )
        self.assertEqual( lsof.select( "AAAAAA", "GGTAAT" ), ("GGTAAT", 2) )

        self.assertEqual( lsof.select( "CCAGTG", "ATCGCA" ), None )
if __name__ == '__main__':
    unittest.main()
