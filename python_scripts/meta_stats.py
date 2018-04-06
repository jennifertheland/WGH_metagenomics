#! /usr/bin/env python


def main():

    import pysam
    from collections import Counter
    global phase_blocks_dict, molecules_per_drop_pos, molecules_per_drop, Counter

    ArgumentParser()
    configureLogging('info')

    with pysam.AlignmentFile(args.rmdup_bam, 'rb') as rmdup_bam:

        phase_blocks_dict = dict()
        molecules_per_drop_pos = dict()
        molecules_per_drop = dict()

        for read in rmdup_bam.fetch('CP010816.1'):

            bc_id = read.get_tag('RG')
            start = read.get_reference_positions()[0]
            stop = read.get_reference_positions()[-1]


            fragments_per_droplet(bc_id,start, stop)

            phase_block_lengths(bc_id, start, stop)


def fragments_per_droplet(bc_id,start,stop):

    try:

        last_stop = molecules_per_drop_pos[bc_id][-1] # all positions in current phase block

        if int(start) - int(last_stop) > 100000:

            try:

                molecules_per_drop[bc_id] += 1

            except KeyError:

                molecules_per_drop[bc_id] = 2

        molecules_per_drop_pos[bc_id] = [start, stop] # reset with new phase block coordinates

    except KeyError:

        molecules_per_drop_pos[bc_id] = [start, stop]


    number_of_molecules_per_drop_list = []

    for key in molecules_per_drop.keys():

        number_of_molecules_per_drop_list.append(molecules_per_drop[key])

    summary_dict_molecules = Counter(number_of_molecules_per_drop_list)

    with open(args.prefix_stats + '.molecules_per_droplet', 'w') as out:

        for key in summary_dict_molecules.keys():

            occurances = summary_dict_molecules[key]

            out.write(str(key) + '\t' + str(occurances) + '\n')


def phase_block_lengths(bc_id, start, stop):
    try:
        start_stop_list = phase_blocks_dict[bc_id]

        if int(start) < int(start_stop_list[0]) and int(start_stop_list[0]) - int(start) < 100000:
            phase_blocks_dict[bc_id][0] = start

        if int(stop) > int(start_stop_list[-1]) and int(stop) - int(start_stop_list[-1]) < 100000:
            phase_blocks_dict[bc_id][-1] = stop

    except KeyError:

        phase_blocks_dict[bc_id] = [start, stop]


    len_of_phase_block_list = []

    for key in phase_blocks_dict.keys():
        start_phase = phase_blocks_dict[key][0]
        stop_phase = phase_blocks_dict[key][-1]

        len_of_phase_block = int(stop_phase) - int(start_phase)

        len_of_phase_block_list.append(len_of_phase_block)

    summary_dict_phase_blocks = Counter(len_of_phase_block_list)

    with open(args.prefix_stats + '.phase_block_lengths', 'w') as out:
        for key in summary_dict_phase_blocks.keys():
            occurances = summary_dict_phase_blocks[key]

            out.write(str(key) + '\t' + str(occurances) + '\n')


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

        parser.add_argument('rmdup_bam',help='Bamfile mapped to reference genome, PCR duplicates removed')

        parser.add_argument('prefix_stats',help='Prefix for the stat files')

#        parser.add_argument('name_3',help='This is an explaination of the argument.')

        args = parser.parse_args()

if __name__ == '__main__': main()