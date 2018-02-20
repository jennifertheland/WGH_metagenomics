#! /usr/bin/env python


def main():

    import sys, time, pysam
    global sys, time

    ArgumentParser()
    configureLogging('info')

    infile = pysam.AlignmentFile(args.name_1, 'rb')

    barcode_id_to_number_of_reads=dict()

    for read in infile.fetch(until_eof=True):

        barcode_id = read.split()[-9]
        try:
            number_of_reads = barcode_id_to_number_of_reads[barcode_id]
            barcode_id_to_number_of_reads[barcode_id] = number_of_reads + 1
        except KeyError:
            barcode_id_to_number_of_reads[barcode_id] = 1

    infile.close()

    infile = pysam.AlignmentFile(args.name_1, 'rb')
    out = pysam.AlignmentFile(args.name_2, 'wb', template=infile)

    for read in infile:

        barcode_id = read.split()[-9]

        if barcode_id_to_number_of_reads[barcode_id] > 20:
            out.write(read)

    infile.close()
    out.close()

def lineCounter(infile):

    counter = 0
    with open(infile,'r') as inf:
        for line in inf:

            counter += 1

    return counter


class Progress():

    def __init__(self,name_of_progress,min,max):

        sys.stderr.write('--- '+ name_of_progress + '---' + '\n')

        self.name_of_progress = name_of_progress
        self.min = min
        self.max = max
        self.one_percent = (self.max-self.min)/100
        self.current_progress = 0
        self.threshold = 0.01

    def progressBarUpdater(self,current_position):

        self.current_progress = current_position/self.max

        if self.current_progress > self.threshold:
            sys.stderr.write('#')
            sys.stderr.flush()
            time.sleep(0.001)
            self.threshold += 0.01

    def terminteProgressbar(self):

        sys.stderr.write('\n')


class configureLogging(object):
    """Configures the logging module. Includes all types of loggind, i.e. INFO, DEBUG,
     WARNING, ERROR and CRITICAL"""

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


class ArgumentParser():

    def __init__(self):

        ArgumentParser.pars(self)

    def pars(self):

        import argparse

        global args

        parser = argparse.ArgumentParser(description=__doc__)

        parser.add_argument('name_1',help='This is an explaination of the argument.')

        parser.add_argument('name_2',help='This is an explaination of the argument.')

        args = parser.parse_args()

if __name__ == '__main__': main()