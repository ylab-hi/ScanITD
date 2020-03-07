import argparse
import sys
import random
import os.path
from StringIO import StringIO

import variation


##
# Generates a donor genome and writes it to the given output path.
#
# @param normal_genome Sequence of the normal genome.
# @param variations List of variations as returned by read_variations.
# @param output_path Donor genome will be written here.
#
def write_donor_genome(normal_genome, variations, output_path):
    mutated_genome = variation.create_indel_genome( normal_genome, variations )
    mutated_genome_file = open( output_path, "w" )
    mutated_genome_file.write( ">mutated-donor\n" )
    mutated_genome_file.write( mutated_genome )
    mutated_genome_file.write( '\n' )
    mutated_genome_file.close( )

##
# Reads the first entry of a fasta file and returns
# its sequence.
#
# @param genome_file A fasta file containing a genome.
#
# @return The sequence of the genome.
#
def read_genome(genome_file):
    # Ignore header
    header = next( genome_file )
    
    genome = StringIO( )
    for line in genome_file:
        if line[0] == '>':
            break

        genome.write( line.strip( ) )

    return genome.getvalue( )

##
# Reads the variations defined in a file and returns
# a them as a list of variation objects.
#
# @param variation_file A file that defines variations.
#
# @return List of variation objects.
#
def read_variations(variation_file):
    variations = [ ]
    for line_number, line in enumerate( variation_file ):
        column = line.split( )
        type = column[ 0 ]
        values = [ int( v ) for v in column[ 1: ] ]

        if column[ 0 ] == "insertion":
            variations.append( variation.Insertion( values[ 0 ], values[ 1 ], -1 ) )
        elif column[ 0 ] == "deletion":
            variations.append( variation.Deletion( values[ 0 ], values[ 1 ] ) )
        elif column[ 0 ] == "duplication":
            variations.append( variation.Insertion( values[ 2 ], values[ 1 ], values[ 0 ] ) )
        elif column[ 0 ] == "translocation":
            variations.append( variation.Deletion( values[ 0 ], values[ 1 ] ) )
            variations.append( variation.Insertion( values[ 2 ], values[ 1 ], values[ 0 ] ) )
        else:
            print( "warning: Ignored bad variation on line: {0}".format( line_number ) )

    return variations

USAGE = """Usage: create_indel_genome genome_file variation_file output_file"""

VARIATION_USAGE = """Path to the variation file. The format of variation_file are lines consisting of the following: 
insertion start length, deletion start length, duplication start length to or translocation start length to."""

if __name__ == '__main__':
    parser = argparse.ArgumentParser( description=USAGE )
    parser.add_argument( 'genome_file', type=argparse.FileType( "r" ), help='Path to the normal genome file.' )
    parser.add_argument( 'variation_file', type=argparse.FileType( "r" ), help=VARIATION_USAGE )
    parser.add_argument( 'output_path', type=str, help='Output path of the mutated genome.' )
    args = parser.parse_args( )
    
    normal_genome = read_genome( args.genome_file )
    variations = read_variations( args.variation_file )
    write_donor_genome( normal_genome, variations, args.output_path ) 
