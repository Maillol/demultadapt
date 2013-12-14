#!/usr/bin/env python
#-*- coding:utf-8 -*-

"""
NAME
       demultadapt - Demultiplex fastq sequences according to tags

SYNOPSIS
       python demulttag.py [OPTION] TAG_FILE

DESCRIPTION
       For each sequence, demultadapt chose an output file and remove the tag.
       
TAG_FILE format
            file with tag sequences i.e ACGTAGCT and output file names, the last line must have * tag :

                tag_1    output_file_name_1
                tag_2    output_file_name_2
                    [...]
                tag_n    output_file_name_3
                *        trash
"""

import sys, os
from davem_fastq import Fastq_read, Fastq_file
import argparse
from itertools import izip
from bisect import bisect_left
 
class FastqFileType( object ) :
    """
    Fastq file factory
    """    
    def __init__( self, mode ) :
        self.mode = mode
    
    def __call__( self, path_name) :
        return Fastq_file( path_name, self.mode )


class Selector( object ) :
    """
    Abstract class to look for an output file in tags_table. 
    When you call select method, this class select an output depending the sequence.
    You must implement __single_select and __paired_select method when subclassing this class.
    """
    def __init__(self, tags_table, single_end) :
        """
        tags_table - see make_tag_table.
        If single_end is True, a monSelector.select( sequence ) call, will call _single_select method
        else a monSelector.select( sequence-1, sequence-2 ) call, will call _paired_select method
        """
        self.tags_table = tags_table
        if single_end :
            self.select = self._single_select
        else :
            self.select = self._paired_select

    def _single_select( self, sequence ) :
        """
        Search line in the tags_table with a sequence
        """
        raise NotImplementedError

    def _paired_select( self, sequence_1, sequence_2 ) :
        """
        Search line in the tags_table with two sequences
        """
        raise NotImplementedError


class Levenshtein_selector( Selector ) :
    tags_table = None
    single_end = False
    rate = 0
    def __init__( self, tags_table, single_end, rate ) :
        if not isinstance( rate, float ) :
            raise ValueError( "rate argument must be a float not %s" % type( rate ) )
        Selector.__init__( self, tags_table, single_end)
        self.rate = rate

    def _single_select( self, sequence) :
        from Levenshtein import ratio
        
        distances = []
        for (adaptator, output_file) in self.tags_table :
            dist = ratio( adaptator, sequence[ : len( adaptator ) ] )
            if dist == 1.0 :
                return (adaptator, output_file)
        
            distances.append( dist )
        
        max_dist = max( distances )
        if max_dist >= self.rate and distances.count( max_dist ) == 1 :
            return self.tags_table[ distances.index( max_dist ) ]
        
        return None

    def _paired_select( self, sequence_1, sequence_2) :
        from Levenshtein import ratio
        distances_1 = []
        distances_2 = []

        for line in self.tags_table :
            adaptator = line[ 0 ]
            dist_1 = ratio( adaptator, sequence_1[ : len( adaptator ) ] )
            dist_2 = ratio( adaptator, sequence_2[ : len( adaptator ) ] )        
            distances_1.append( dist_1 )
            distances_2.append( dist_2 )

        max_dist_1 = max( distances_1 )
        max_dist_2 = max( distances_2 )

        if max_dist_1 > max_dist_2 :
            if max_dist_1 >= self.rate and distances_1.count( max_dist_1 ) == 1 :
                return self.tags_table[ distances_1.index( max_dist_1 ) ]
                                
        elif max_dist_1 < max_dist_2 :
            if max_dist_2 >= self.rate and distances_2.count( max_dist_2 ) == 1 :
                return self.tags_table[ distances_2.index( max_dist_2 ) ]

        else :
            if max_dist_1 >= self.rate :
                if distances_1.count( max_dist_1 ) == 1 :
                    index_1 = distances_1.index( max_dist_1 )
                    index_2 = distances_2.index( max_dist_2 )
                    if index_1 == index_2 :
                        return self.tags_table[ index_1 ]

                elif distances_2.count( max_dist_2 ) == 1 :
                    index_1 = distances_1.index( max_dist_1 )
                    index_2 = distances_2.index( max_dist_2 )
                    if index_1 == index_2 :
                        return self.tags_table[ distances_2.index( max_dist_2 ) ]

        return None


class LevenshteinAllSelector( Levenshtein_selector ) :
    """
    Levenshtein_selector with paired-end way, Two member sequence of paired-end must be greater or equal to 
    rate min and have same tags.
    """

    def _paired_select( self, sequence_1, sequence_2) :
        from Levenshtein import ratio
        distances_1 = []
        distances_2 = []

        for line in self.tags_table :
            adaptator = line[ 0 ]
            dist_1 = ratio( adaptator, sequence_1[ : len( adaptator ) ] )
            dist_2 = ratio( adaptator, sequence_2[ : len( adaptator ) ] )        
            distances_1.append( dist_1 )
            distances_2.append( dist_2 )

        max_dist_1 = max( distances_1 )
        max_dist_2 = max( distances_2 )

        if ( max_dist_1 >= self.rate and max_dist_2 >= self.rate 
           and distances_1.count( max_dist_1 ) == distances_2.count( max_dist_2 ) == 1 ) :
               adapt_1 = self.tags_table[ distances_1.index( max_dist_1 ) ] 
               adapt_2 = self.tags_table[ distances_2.index( max_dist_2 ) ]
               if adapt_1 == adapt_2 :
                    return adapt_1
        else :
            return None

class Std_selector( Selector ):
    """
    Search in the tags_table, sequence start must be identical to tag not similar. 
    """

    def _paired_select( self, sequence_1, sequence_2):
        l1 = self._single_select( sequence_1 )
        l2 = self._single_select( sequence_2 )
        if l1 is None :
            return l2
        
        if l2 is None :
            return l1
        
        if l1 == l2 :
            return l1

        return None
        
        
    def _single_select( self, sequence):
        a = 0
        b = len( self.tags_table ) -1
        if b == -1 :
            return None
        
        while a <= b  :
            m = ( a + b ) // 2
            adaptator = self.tags_table[ m ][ 0 ]
            start_seq = sequence[ : len( adaptator ) ]
            
            if adaptator > start_seq :
                b = m - 1
            elif adaptator < start_seq :
                a = m + 1
            else :
                return self.tags_table[ m ]

        if adaptator == sequence[ : len( adaptator ) ] :
            return self.tags_table[ m ]
        return None


def get_adapt_counter( opened_adapt_file ) :
    """
    return { tag1 : 0, 
             tag2 : 0,
             ... 
             tagN : 0 }
    """
    d = {}
    opened_adapt_file.seek(0)
    for line in opened_adapt_file :
        if not line.isspace() :
            try :    
                adapt, name_tag = line.split()
            except ValueError :
                print >> sys.stderr, "Error '%s' is an invalid file format." %  opened_adapt_file.name
                exit( 1 )
            d[ adapt ] = [ name_tag, 0 ]
    return d


def get_maximal_annalogie( file_adapt ) :
    """
    Compute the max Levenshtein rate between tags.
    """
    from Levenshtein import ratio
    adaptators = []
    for line in file_adapt :
        if line :
            (adapt, name) = line.split()
            if adapt != "*" :           
                adaptators.append( adapt )
   
    ratio_max = 0.0
    for i, adapt in enumerate( adaptators ) :
        for adapt2 in adaptators[i+1:] :
            ratio_max = max( ratio_max,ratio( adapt, adapt2 ) )

    return ratio_max




def make_tag_table( opened_adapt_file, prefix, paired_end=True ) :
    """
    Return the output file list (tag_table) and trash_file. 
    
    When paired_end is True, the tag_table format is:
    
        [ (tagA, output_file_A_1, output_file_A_2 ),
          (tagB, output_file_B_1, output_file_B_2 ),
          ...  ]
    
        return ( tag_table, (trash_file_1, trash_file_2) )
    
    else:
        [ (tagA, output_file_A ),
          (tagB, output_file_B ), 
          ...  ]
        return ( tag_table, (trash_file_1,) )
    """ 
    

    ada_files = []
    default = None
    cache_name_file_by_adapt = {}

    for line in opened_adapt_file :
        if not line.isspace() :
                try :    
                    adapt, suffix_file = line.split()
                except ValueError :
                    print >> sys.stderr, "Error: '%s' is an invalid file format." %  opened_adapt_file.name
                    exit( 1 )

                if paired_end :
                    if line[0] == '*' :
                        default = ( Fastq_file( "%s-%s_1.fastq" % (prefix, suffix_file), "w" ),
                                    Fastq_file( "%s-%s_2.fastq" % (prefix, suffix_file), "w" ), )

                    else :
                        if suffix_file in cache_name_file_by_adapt :
                            f1, f2 = cache_name_file_by_adapt[ suffix_file ]
                            ada_files.append( ( adapt, f1, f2 ) )

                        else :
                            f1 = Fastq_file( "%s-%s_1.fastq" % (prefix, suffix_file), "w" )
                            f2 = Fastq_file( "%s-%s_2.fastq" % (prefix, suffix_file), "w" )
                            ada_files.append( (adapt, f1, f2) )
                            cache_name_file_by_adapt[ suffix_file ] = (f1, f2)


                else :
                    # TODO faire le system de cache pour le mode single.
                    if line[0] == '*' :
                        default = ( Fastq_file( "%s-%s.fastq" % (prefix, suffix_file), "w" ) , )

                    else :
                        ada_files.append(  (
                                                adapt, 
                                                Fastq_file( "%s-%s.fastq" % (prefix, suffix_file), "w" )
                                           )
                        )

    if default is None :
        print >> sys.stderr, "Le fichier '%s' n'a pas de ligne avec le tag jocker *.\nAjouter une ligne '*    tag_name'." %  opened_adapt_file.name
        sys.exit(1)
    
    ada_files.sort()
    return ada_files, default


def parse_user_argument() :
    """
    Get user argument.
    """
    parser = argparse.ArgumentParser( formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__ )
    parser.add_argument( 'file_adapt', metavar="FILE_TAG", nargs=1, type=argparse.FileType('r') )

    parser.add_argument( '-f', '--fastq_1', dest="fastq_1", type=FastqFileType( "r" ), action='store', 
                            help="single-end file or paired-end file 1" )

    parser.add_argument( '-F', '--fastq_2', dest="fastq_2", type=FastqFileType( "r" ), action='store', default=None,
                            help="paired-end file 2" )

    parser.add_argument( '-p', '--output_prefix', dest="output_prefix", default="", action='store',
                            help="output file names are: PREFIX-NAME_IN_FILE_TAG.fastq"  )

    parser.add_argument( '-l', '--levenshtein', dest="levenshtein", action='store', type=float, default=None,
                            help="Use a Levenshtein distance to demultiple" )

    parser.add_argument( '-v', '--verbose', dest="verbose", action='store_true',
                            help="explain what is being done" )

    parser.add_argument( '-a', '--analogy', dest="analogy", action='store_true',
                            help="Compute the maximal Levenshtein ratio between adaptors" )

    parser.add_argument( '--all', dest="all", action='store_true',
                            help="if is enable, and levenshtein too with paired-end mode, All members of the paired-end must have rate greater than or equal to levenshtein rate and the same tag." )

    user_args = parser.parse_args()
    user_args.file_adapt = user_args.file_adapt[0]
    user_args.single_end = user_args.fastq_2 is None 
    return user_args

def main() :
    user_args = parse_user_argument()

    if user_args.analogy :
        print "Maximal Levenshtein ratio between adaptors is %f" % get_maximal_annalogie( user_args.file_adapt )
        sys.exit(0)        
            
    output_files_by_adapt, defaults_files = make_tag_table( user_args.file_adapt,
                                                            user_args.output_prefix,
                                                            not user_args.single_end )

    nb_reads_writen = get_adapt_counter( user_args.file_adapt )

    user_args.file_adapt.close()

    if user_args.levenshtein : 
        if user_args.all :
            select_output_file = LevenshteinAllSelector( output_files_by_adapt, 
                                                       user_args.single_end,
                                                       user_args.levenshtein )

        else :
            select_output_file = Levenshtein_selector( output_files_by_adapt, 
                                                       user_args.single_end,
                                                       user_args.levenshtein )
    else :
        select_output_file = Std_selector( output_files_by_adapt, 
                                           user_args.single_end )

    # Single_end
    if user_args.single_end :
        default_file = defaults_files[0]
        for str_read in user_args.fastq_1 :
            read = Fastq_read( str_read )
            adapt_and_line = select_output_file.select( read.seq )
            if adapt_and_line is None :
                if user_args.verbose :
                    print "Read '%s' start with %s... and go to *" % (read.name, read.seq[ : 14 ])
                default_file.write( str( read ) )
                nb_reads_writen[ '*' ][ 1 ] += 1

            else :
                (adapt, output_file) = adapt_and_line
                if user_args.verbose :
                    print "Read '%s' start with %s... and go to %s" % (read.name, read.seq[ : len( adapt ) ], adapt)

                read.cut_start( len( adapt ) )
                output_file.write( str( read ) )
                nb_reads_writen[ adapt ][ 1 ] += 1

        user_args.fastq_1.close()

        for adapt, output_file in output_files_by_adapt :
            output_file.close()

    # Paired_end
    else :
        (default_file_1, default_file_2) = defaults_files

        for str_read_1, str_read_2 in izip( user_args.fastq_1, user_args.fastq_2 ) :
            read_1 = Fastq_read( str_read_1 )
            read_2 = Fastq_read( str_read_2 )

            adapt_and_line = select_output_file.select( read_1.seq, read_2.seq )
            
            if adapt_and_line is None :
                default_file_1.write( str( read_1 ) )
                default_file_2.write( str( read_2 ) )
                nb_reads_writen[ '*' ][1] += 1

            else :
                (adapt, output_file_1, output_file_2 ) = adapt_and_line

                read_1.cut_start( len( adapt ) )
                read_2.cut_start( len( adapt ) )

                output_file_1.write( str( read_1 ) )
                output_file_2.write( str( read_2 ) )
                nb_reads_writen[ adapt ][1] += 1

        user_args.fastq_1.close()
        user_args.fastq_2.close()

        for adapt, file_1, file_2 in output_files_by_adapt :
            file_1.close()
            file_2.close()

    # show stat.
    for nb_reads_by_name in nb_reads_writen.values() :
        print "%s %d reads" % tuple( nb_reads_by_name )


if __name__ == '__main__':
    main()
