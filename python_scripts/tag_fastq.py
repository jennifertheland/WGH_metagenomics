#! /usr/bin/env python

def main():
    """ Takes a .fastq file with trimmed reads and a .fasta file with clustered barcodes. Output is a tagged .fastq file. """

    global args, summaryInstance
    from Bio import SeqIO
    import logging



    configureLogging('info')

    readArgs() # Argument parsing
    logging.info('Arguments read successfully')

    summaryInstance = Summary()

    # Generates a dictionary from concatenated .clstr file (stores in Summary instance)
    # dict1 = consensus barcode representing the cluster :  non-consensus barcode in the same cluster
    # dict2 = master dictionary, concatenated from dict1

    readAndProcessClusters(args.input_clstr)

    logging.info('Cluster file processed successfully')

    # Tagging bam mapping entries with RG:Z:bc_seq BC:Z:bc_seq
    with open(args.input_trimmed_fastq, 'r') as openin, open(args.output_tagged_fastq, 'w') as openout:

        for record in SeqIO.parse(openin,'fastq'):

            bc_seq = record.id.split('_')[-1]

            try:
                consensus_seq = summaryInstance.master_barcode_dict[bc_seq]
            except KeyError:
                continue

            record.description = record.id + ' BC:Z:' + consensus_seq + '-1'

            SeqIO.write(record,openout,'fastq')

    logging.info('Tagging completed')

def readAndProcessClusters(openClstrFile):
    """ Reads clstr file and builds consensus_bc:bc dict in Summary instance."""
    new_cluster = False

    with open(openClstrFile,'r') as f:

        clusterInstance = ClusterObject(clusterId=f.readline()) # set clusterId for the first cluster

        for line in f:
            # Reports cluster to master dict and start new cluster instance

            if line.startswith('>'):
                new_cluster = True

            elif new_cluster == True:
                new_cluster = False
                clusterInstance = ClusterObject(clusterId=line)
                clusterInstance.addBarcodeToDict(line)
                summaryInstance.updateReadToClusterDict(clusterInstance)

            else:
                clusterInstance.addBarcodeToDict(line)

        # Add last cluster to master dict
        summaryInstance.updateReadToClusterDict(clusterInstance)

class ClusterObject(object):
    """ Cluster object, one for each cluster"""

    def __init__(self, clusterId):

        self.consensus_to_bc_dict = dict()
        self.consensus = clusterId.split()[-2].split(':')[-1].rstrip('.')

    def addBarcodeToDict(self, line):
        '''Add non-consensus bc to consensus_to_bc_dict'''
        accession = line.split()[2].rstrip('.')
        barcode = accession.split(':')[-1]
        self.consensus_to_bc_dict[barcode] = self.consensus

class readArgs(object):
    """ Reads arguments and handles basic error handling like python version control etc."""

    def __init__(self):
        """ Main funcion for overview of what is run. """

        readArgs.parse(self)
        readArgs.pythonVersion(self)

    def parse(self):

        import argparse, multiprocessing
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("input_trimmed_fastq", help=".fastq file with trimmed reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
        parser.add_argument("output_tagged_fastq", help=".fastq file tagged with barcode consensus sequence")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. Example: python "
                                 "TagGD_prep.py -p 2 insert_r1.fq unique.fa")

        args = parser.parse_args()

    def pythonVersion(self):
        """ Makes sure the user is running python 3."""

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

    def processors(self):

        import multiprocessing
        processor_count = args.processors
        max_processor_count = multiprocessing.cpu_count()
        if processor_count == max_processor_count:
            pass
        elif processor_count > max_processor_count:
            sys.stderr.write(
                'Computer does not have ' + str(processor_count) + ' processors, running with default (' + str(
                    max_processor_count) + ')\n')
            processor_count = max_processor_count
        else:
            sys.stderr.write('Running with ' + str(processor_count) + ' processors.\n')

        return processor_count

class Summary(object):
    """ Summarizes chunks"""

    def __init__(self):
        self.master_barcode_dict = dict()

    def updateReadToClusterDict(self, input_object):
        """ Merges cluster-specific dictionaries to a master dictionary."""

        input_dict = input_object.consensus_to_bc_dict

        for barcode in input_dict.keys():
            self.master_barcode_dict[barcode] = input_object.consensus

class configureLogging(object):

    import logging
    global logging

    def __init__(self,logType):

        log_format = ' %(asctime)s %(levelname)s: < %(message)s >'

        if logType == 'info':
            configureLogging.INFO(self,log_format)
        elif logType == 'warning':
            configureLogging.WARNING(self,log_format)
        elif logType == 'error':
            configureLogging.ERROR(self,log_format)
        elif logType == 'critical':
            configureLogging.CRITICAL(self,log_format)
        elif logType == 'debug':
            configureLogging.DEBUG(self,log_format)

    def INFO(self,log_format):

        logging.basicConfig(level=logging.INFO, format=log_format)

    def WARNING(self,log_format):

        logging.basicConfig(level=logging.WARNING, format=log_format)

    def ERROR(self,log_format):

        logging.basicConfig(level=logging.ERROR, format=log_format)

    def CRITICAL(self,log_format):

        logging.basicConfig(level=logging.CRITICAL, format=log_format)

    def DEBUG(self,log_format):

        logging.basicConfig(level=logging.DEBUG, format=log_format)


if __name__=="__main__": main()
