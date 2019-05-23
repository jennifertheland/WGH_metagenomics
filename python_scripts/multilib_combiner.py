#! /usr/bin/env python2

def main():

    #
    # Imports & globals
    #
    global args, summaryInstance, sys, time
    import pysam, sys, time, BLR_functions as BLR

    #
    # Argument parsing
    #
    argumentsInstance = readArgs()

    #
    # Process data
    #

    # Error search: If not equally dividable by 3 => exit
    if len(args.read_files) & 1:
        if not args.force_run:
            sys.exit('ARGUMENT ERROR. Not an even number of arguments, please supply both R1 and R2 for all read groups')

    BLR.report_progress("Starting analysis")

    # Initials
    bc_set = int()
    with open(args.r1_out, 'w') as r1_out, open(args.r2_out, 'w') as r2_out:
        read_files = list()
        for read_file in args.read_files:

            # Step two files at the time
            read_files.append(read_file)
            if not len(read_files) == 2:
                continue
            else:
                bc_set += 1
                BLR.report_progress('Starting with RG:\t' + str(bc_set))

            # Setup variables
            read1_file = read_files[0]
            read2_file = read_files[1]
            progress = BLR.ProgressReporter('Reading files:\t' + str(read1_file) + '\t' + str(read2_file), 1000000)
            generator = BLR.FileReader(read1_file, read2_file)
            for read1, read2 in generator.fastqPairedReader():

                # Fetch bc seq
                header = read1.header.split()

                # Add read group to header: @HEADER_bcSeq RG:Z:rg-N BC:Z:bc_clstr_id-N (where N is rg identifier)
                read1.header = header[0] + '\tRG:Z:rg-' + str(bc_set) + '\t' + header[1][:-2] + '-' + str(bc_set)
                read2.header = read1.header

                # Write to out
                r1_out.write(read1.fastq_string())
                r2_out.write(read2.fastq_string())
                progress.update()

            # Empty variables
            generator.close()
            read_files = list()
    BLR.report_progress("Finished")

class readArgs(object):
    """
    Reads arguments and handles basic error handling like python version control etc.
    """

    def __init__(self):

        readArgs.parse(self)
        readArgs.pythonVersion(self)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("read_files", nargs='+' ,help="read files")
        parser.add_argument("-r1", "--r1_out", required=True ,help="Read one out")
        parser.add_argument("-r2", "--r2_out", required=True, help="Read two out")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-s", "--split", type=str, default='_', help="Character splitting the header from the "
                                                                         "barcode sequence and cluster ID.")

        args = parser.parse_args()

    def pythonVersion(self):
        """ Makes sure the user is running python 3."""

        #
        # Version control
        #
        import sys
        if sys.version_info.major == 3:
            pass
        else:
            sys.stderr.write('\nWARNING: you are running python ' + str(
                sys.version_info.major) + ', this script is written for python 3.')
            if not args.force_run:
                sys.stderr.write('\nAborting analysis. Use -F (--Force) to run anyway.\n')
                sys.exit()
            else:
                sys.stderr.write('\nForcing run. This might yield inaccurate results.\n')

if __name__=="__main__": main()