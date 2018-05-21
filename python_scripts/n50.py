#!/usr/bin/env python3


def main():

    from Bio import SeqIO
    global args

    ArgumentParser()

    with open(args.name_1, 'r') as openin:

        len_of_all = []

        for record in SeqIO.parse(openin,'fasta'):

            seq_len_temp = len(record.seq)

            len_of_all.append(seq_len_temp)

        total_length = sum(len_of_all)
        len_of_all_sort = sorted(len_of_all)
        cumulative_len = 0
        for i in len_of_all_sort:

            cumulative_len += i

            if cumulative_len/total_length > 0.5:

                print('N50: %d' % i)

                break



class ArgumentParser():

    def __init__(self):

        ArgumentParser.parse(self)

    def parse(self):

        import argparse

        global args

        parser = argparse.ArgumentParser(description=__doc__)

        parser.add_argument('name_1',help='This is an explaination of the argument.')

        args = parser.parse_args()

if __name__ == '__main__': main()
