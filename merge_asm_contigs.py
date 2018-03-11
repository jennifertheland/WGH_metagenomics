#! /usr/bin/env python

def main():

    import pysam
    global args

    ArgumentParser()
    configureLogging('info')

    in_bam = pysam.AlignmentFile(args.final_asm_mapped, "rb",check_sq=False)

    bc_to_contig = dict()

    for read in in_bam.fetch(until_eof=True):

        try:
            if read.reference_id not in bc_to_contig[read.get_tag('BC')] and read.reference_id != -1:

                bc_to_contig[read.get_tag('BC')].append(read.reference_id)

        except KeyError:

            if read.reference_id != -1:
                bc_to_contig[read.get_tag('BC')] = [read.reference_id]

    logging.info('Dictionary created successfully')

    in_bam.close()

    with open(args.output_asm,'w') as out:

        spanning_bcs = dict()

        for key in bc_to_contig.keys():

            if len(bc_to_contig[key]) > 1:

                spanning_bcs[key] = bc_to_contig[key]

    print(spanning_bcs)

class ArgumentParser():

    def __init__(self):

        ArgumentParser.pars(self)

    def pars(self):

        import argparse

        global args

        parser = argparse.ArgumentParser(description=__doc__)

        parser.add_argument('final_asm_mapped',help='This is an explaination of the argument.')

        parser.add_argument('output_asm', help='This is an explaination of the argument.')

        args = parser.parse_args()

class configureLogging(object):
    """Configures the logging module. Includes all types of loggind, i.e. INFO, DEBUG,
     WARNING, ERROR and CRITICAL"""

    import logging
    global logging

    def __init__(self, logType):

        log_format = ' %(asctime)s %(levelname)s: < %(message)s >'

        if logType == 'info':
            configureLogging.INFO(self, log_format)
        elif logType == 'warning':
            configureLogging.WARNING(self, log_format)
        elif logType == 'error':
            configureLogging.ERROR(self, log_format)
        elif logType == 'critical':
            configureLogging.CRITICAL(self, log_format)
        elif logType == 'debug':
            configureLogging.DEBUG(self, log_format)

    def INFO(self, log_format):

        logging.basicConfig(level=logging.INFO, format=log_format)

    def WARNING(self, log_format):

        logging.basicConfig(level=logging.WARNING, format=log_format)

    def ERROR(self, log_format):

        logging.basicConfig(level=logging.ERROR, format=log_format)

    def CRITICAL(self, log_format):

        logging.basicConfig(level=logging.CRITICAL, format=log_format)

    def DEBUG(self, log_format):

        logging.basicConfig(level=logging.DEBUG, format=log_format)

if __name__ == '__main__': main()