#!/usr/bin/env python
#-*- coding:utf-8 -*-
import sys, os
sys.path.append( "../" )

from davem_fastq import Fastq_file, Fastq_read

def check_returned_single_file() :
    for file_name in os.listdir("./returned") :
        fq = Fastq_file( "./returned/" + file_name, "r" )
        for read in fq :
            if Fastq_read( read ).name[6:] != file_name[2:-6]:
                print >> sys.stderr, "Error, read %s is in file %s" % (Fastq_read( read ).name, file_name)
                sys.exit( 1 )
        
        print "Contain of file %-20s      [OK]" % file_name
        os.remove( "./returned/" + file_name )



def check_returned_paired_end_file() :
    for file_name in os.listdir("./returned") :
        fq = Fastq_file( "./returned/" + file_name, "r" )
        for read in fq :
            print >> sys.stderr, "$$$", Fastq_read( read ).name, file_name
            if Fastq_read( read ).name[8:] != file_name[2:-8]:
                print >> sys.stderr, "Error, read %s is in file %s" % (Fastq_read( read ).name, file_name)
                sys.exit( 1 )
        
        print "Contain of file %-20s      [OK]" % file_name
        os.remove( "./returned/" + file_name )



print "Test single levenshtein 1"
os.system( "python ../demultadapt.py -v -f single-l1.fastq -p returned/r -l 0.80 adapt.txt" )
check_returned_single_file()

print "Test single levenshtein 2"
os.system( "python ../demultadapt.py  -f single-strcmp.fastq -p returned/r -l 1.0 adapt.txt" )
check_returned_single_file()


print "Test single 100%"
os.system( "python ../demultadapt.py  -f single-strcmp.fastq -p returned/r  adapt.txt" )
check_returned_single_file()

print "Test paired levenshtein"
os.system( "python ../demultadapt.py  -f paired-l1-1.fastq -F paired-l1-2.fastq -l 0.80 -p returned/r  adapt.txt" )
check_returned_paired_end_file()


print "All test pass"
