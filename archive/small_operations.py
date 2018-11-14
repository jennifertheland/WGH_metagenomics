#! /usr/bin/env python


def main():

    ArgumentParser()
    configureLogging('info')

    with open(args.name_1,'r') as inf, open(args.name_2,'w') as out:

        #step_size = 5000

        number_of_bins = 1

        for line in inf:

            if number_of_bins == 1:

                name_chr1, start_pos1, stop_pos1, coverage1 = line.split()

                number_of_bins += 1

            #if int(phase_block_length) < step_size * number_of_bins and int(phase_block_length) > step_size * number_of_bins - 5000:

            #   total_occurances += int(occurance)

            elif number_of_bins == 2:

                name_chr2, start2_pos, stop_pos2, coverage2 = line.split()

                total_cov = (float(coverage1) + float(coverage2))/2

                out.write(name_chr1 + '\t' + start_pos1 + '\t' + stop_pos2 + '\t' + str(total_cov) + '\n')

                number_of_bins = 1

            #    out.write(str(step_size * number_of_bins - 5000) + '\t' + str(step_size*number_of_bins) + '\t' + str(total_occurances) + '\n')

            #   if int(phase_block_length) > step_size * number_of_bins:

            #        while int(phase_block_length) > step_size * number_of_bins:


             #           number_of_bins += 1

              #      total_occurances = int(occurance)

             #   else:
              #      number_of_bins += 1

        #out.write(str(step_size * number_of_bins - 5000) + '\t' + str(step_size*number_of_bins) + '\t' + str(total_occurances) + '\n')


        # counter = 0
        # total_coverage = 0
        #
        # for line in inf:
        #
        #     if counter == 0:
        #
        #         first_position = line.split()[1]
        #
        #     counter += 1
        #
        #     pos_coverage = line.split()[-1]
        #
        #     total_coverage += int(pos_coverage)
        #
        #     if counter == 50:,
        #
        #         col_1,col_2,col_3,col_4 = line.split()
        #
        #         number_of_bases = int(col_3) - int(first_position)
        #
        #         coverage = total_coverage / number_of_bases
        #
        #         out.write('chr1' + '\t' + first_position + '\t' + col_2 + '\t' + str(coverage) + '\n')
        #
        #         counter = 0
        #
        #         total_coverage = 0
        #
        #         coverage = 0

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

#        parser.add_argument('name_3',help='This is an explaination of the argument.')

        args = parser.parse_args()

if __name__ == '__main__': main()