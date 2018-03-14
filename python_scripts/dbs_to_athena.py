#! /usr/bin/env python


def main():

    from Bio import SeqIO

    ArgumentParser()
    configureLogging('info')

    with open(args.name_1,'r') as inf:

        logging.info('File read.')

        id_to_description_dict = dict()

        for record in SeqIO.parse(inf, 'fastq'):

            id_only = record.id.split('_')[0]

            id_to_description_dict[id_only] = record.description

        logging.info('Dictionary created.')

    with open(args.name_2,'r') as inf2, open(args.name_3,'w') as out:

        for record in SeqIO.parse(inf2,'fastq'):

            try:
                full_id = id_to_description_dict[record.id]

                record.id = full_id

                SeqIO.write(record,out,'fastq')

            except KeyError:
                continue

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

        parser.add_argument('name_3', help='This is an explaination of the argument.')

        args = parser.parse_args()

if __name__ == '__main__': main()