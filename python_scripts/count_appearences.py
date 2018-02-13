#! /usr/bin/env python


def main():

    readArgs()

    with open(args.input_txt,'r') as inf, open(args.output_summary,'w') as out:

        previous_line = inf.readline()

        counter = 1

        for line in inf:
            if line == previous_line:
                counter += 1
                previous_line = line
            else:
                out.write(previous_line + '\t' + str(counter) + '\n')
                previous_line =  line
                counter = 1

class readArgs(object):
    """ Reads arguments and handles basic error handling like python version control etc."""

    def __init__(self):
        """ Main funcion for overview of what is run. """

        readArgs.parse(self)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("input_txt", help=".bam file with mapped reads which is to be tagged with barcode id:s.")
        parser.add_argument("output_summary", help=".clstr file from cdhit clustering.")


        args = parser.parse_args()

if __name__ == '__main__': main()