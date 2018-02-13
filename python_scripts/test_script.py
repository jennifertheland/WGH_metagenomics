#! /usr/bin/env python



def main():

    global args

    args_ins = readArgs()

    with open(args.input_rm_file ,'r') as new_id_uniq, open(args.output_bc_id_consensus,'w') as out:

        number_of_processed_lines=0
        id_to_unique_dict, id_to_consensus_dict = createIDConsensusDict(args.input_clstr, args.input_bc_id_unique_dict)
        out.write('new ID \t new bc_seq \t\t old ID \t old bc_seq \t old unique sequence \n')

        for line in new_id_uniq:

            number_of_processed_lines += 1

            new_bc_id, unique_bc_seq = line.split()

            old_bc_id = id_to_unique_dict[unique_bc_seq] # consensus barcode id

            new_bc_seq = id_to_consensus_dict[int(new_bc_id)]

            old_consensus_seq = id_to_consensus_dict[old_bc_id]

            out.write(new_bc_id + '\t' + new_bc_seq + '\t' + str(old_bc_id) + '\t' + old_consensus_seq + '\t' + unique_bc_seq + '\n')

            if number_of_processed_lines % 100 == 0:
                print('Number of processed lines' + str(number_of_processed_lines))


def createIDConsensusDict(infile1,infile2):

    id_to_consensus_dict = dict()
    import ast

    with open(infile2, 'r') as inf2:

        for first_line in inf2:

            id_to_unique_dict = first_line # unique barcodes to cluster id

            id_to_unique_dict = ast.literal_eval(id_to_unique_dict)

    with open(infile1,'r') as inf1:

        for line in inf1:

            if line.startswith('>'):
                continue

            elif line.startswith('0'): # grabs consensus seqs

                consensus_accession = line.split()[2].rstrip('.')  # inget att bry sig om
                consensus_barcode = consensus_accession.split(':')[-1]  # consensus seq för ett kluster
                bc_id = id_to_unique_dict[consensus_barcode]

                id_to_consensus_dict[bc_id] = consensus_barcode
    return id_to_unique_dict, id_to_consensus_dict

class readArgs(object):
    """ Reads arguments and handles basic error handling like python version control etc."""

    def __init__(self):
        """ Main funcion for overview of what is run. """

        readArgs.parse(self)

    def parse(self):

        #
        # Imports & globals
        #
        import argparse, multiprocessing
        global args

        parser = argparse.ArgumentParser(description=__doc__)

        # Arguments
        parser.add_argument("input_bc_id_unique_dict", help=".bam file with mapped reads which is to be tagged with barcode id:s.")
        parser.add_argument("input_rm_file", help=".clstr file from cdhit clustering.")
        parser.add_argument("input_clstr", help=".clstr file from cdhit clustering.")
        parser.add_argument("output_bc_id_consensus", help=".bam file with barcode cluster id in the bc tag.")

        # Options
        parser.add_argument("-F", "--force_run", action="store_true", help="Run analysis even if not running python 3. "
                                                                           "Not recommended due to different function "
                                                                           "names in python 2 and 3.")
        parser.add_argument("-p", "--processors", type=int, default=multiprocessing.cpu_count(),
                            help="Thread analysis in p number of processors. Example: python "
                                 "TagGD_prep.py -p 2 insert_r1.fq unique.fa")
        parser.add_argument("-e", "--exclude_N", type=bool, default=True, help="If True (default), excludes .bam file "
                                                                               "reads with barcodes containing N.")

        args = parser.parse_args()

if __name__ == '__main__': main()