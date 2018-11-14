#! /usr/bin/env python

def main():
    """ Takes a fastq file and sorts it according to barcode sequence """

    global args

    argumentsInstance = readArgs()

    with open(args.one_line_file_records, 'r') as openin, open(args.reformatted_fastq_sorted, 'w') as openout:
        for line in openin:
            arguments = line.split()
            openout.write(arguments[0] + ' ' + arguments[1] + '\n')
            openout.write(arguments[2] + '\n')
            openout.write(arguments[3] + '\n')
            openout.write(arguments[4] + '\n')

class readArgs(object):
    """ Reads arguments and handles basic error handling like python version control etc."""

    def __init__(self):
        """ Main funcion for overview of what is run. """

        readArgs.parse(self)

    def parse(self):

        import argparse
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("one_line_file_records", help="temp file where a .fastq record is in a single line")
        parser.add_argument("reformatted_fastq_sorted", help="outfile, sorted .fastq based on barcode sequence")

        args = parser.parse_args()

if __name__ == '__main__': main()