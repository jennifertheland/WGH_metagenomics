#! /usr/bin/env python


def main():

    import sys, time
    global sys, time

    ArgumentParser()
    configureLogging('info')


    with open(args.name_1,'r') as inf, open(args.name_2,'w') as outf:

        logging.info('files read successfully')

        # number_of_lines = lineCounter(args.name_of_argument)
        # counter = 0
        # pg_bar_read_file = Progress("Reading file ", 1, number_of_lines)

        out_str = inf.readline()
        counter = 0
        sum_of_contigs = 0
        number_of_contigs = 0

        for line in inf:

            if line.startswith('>'):

                number_of_contigs += 1

                if counter != 0:

                    out_str += str(counter) + '\n'
                    outf.write(out_str)
                    sum_of_contigs += counter
                    out_str = ''
                    counter = 0
                    out_str += line

                else:
                    continue

            else:

                counter += len(line)
        mean_len_of_contigs = sum_of_contigs /number_of_contigs
        outf.write('Total length of contigs: ' +str(sum_of_contigs))
        outf.write('Mean length of contigs: ' + str(mean_len_of_contigs) + '\n')

        #     counter += 1
        #     pg_bar_read_file.progressBarUpdater(counter)
        #
        # pg_bar_read_file.terminteProgressbar()

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