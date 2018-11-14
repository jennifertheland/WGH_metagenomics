#! /usr/bin/env python


def main():

    ArgumentParser()
    configureLogging('info')


    with open(args.name_1,'r') as posfile, open(args.name_2,'r') as reference, open(args.name_3,'w') as out:

        """posfile = chromosome name newline pos1 pos2, reference = reference file without linebreaks"""

        dict1 = {}

        for line in reference:

            if line.startswith(('>AM260479.1_chr1','>AM260480.1_chr2','>AY305378.1_megaplasmid','>CP010816.1_BL21_TaKaRa')):

                temp = line.split('>')[-1]
                genome = temp.strip('\n')

            else:

                dict1[genome] = line

        for line in posfile:

            if line.startswith(('AM260479.1_chr1', 'AM260480.1_chr2', 'AY305378.1_megaplasmid', 'CP010816.1_BL21_TaKaRa')):

                genome = line.strip('\n')
                out.write(genome + '\n')

            else:

                pos1, pos2 = line.split()

                sequence = dict1[genome][int(pos1)-1:int(pos2)-1]

                out_str = content_calculator(sequence)

                if len(out_str) != 0:

                    out.write(out_str)






def content_calculator(sequence):

    content = [0,0,0,0] # A, T, G, C

    out_str = ''

    for i in sequence:

        if i == 'A':

            content[0] += 1

        elif i == 'T':

            content[1] += 1

        elif i == 'G':

            content[2] += 1

        elif i == 'C':

            content[3] += 1

    if len(sequence) != 0:

        percentage_content = ["%.2f" % (i/len(sequence)) for i in content]

        gc_content = float(percentage_content[-1]) + float(percentage_content[-2])

        out_str = str(percentage_content) + '\t' + str(gc_content) + '\t' + str(len(sequence)) + '\n'

    return out_str









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

        parser.add_argument('name_1',help='Position file')

        parser.add_argument('name_2',help='Reference')

        parser.add_argument('name_3',help='out')

        args = parser.parse_args()

if __name__ == '__main__': main()