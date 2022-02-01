#!/usr/bin/env python
""" convert alignment information from blast table to dotplot like density matrix"""
import numpy as np
import csv
import argparse
import subprocess
import tempfile


def blast2density(w, length_of_seq, input_blast):
    """
    :param w: width of the array
    :param length_of_seq: length of DNA sequence
    :param input_blast: file with blast output
    :return: density matrix
    """
    mat_ff = np.zeros([w + 1, w + 1], dtype=float)
    mat_rc = np.zeros([w + 1, w + 1], dtype=float)
    i = 0
    with open(input_blast, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t', )
        for line in csvreader:
            i += 1
            if line[0][0] == "#":
                continue
            pid = float(line[2])
            # coordinates of alignment
            x0 = float(line[6])
            x1 = float(line[7])
            y0 = float(line[8])
            y1 = float(line[9])
            # coordinates for matrix
            mx0 = w * x0 / length_of_seq
            mx1 = w * x1 / length_of_seq
            my0 = w * y0 / length_of_seq
            my1 = w * y1 / length_of_seq
            lx = round(abs(mx1 - mx0) + 0.5)
            ly = round(abs(my1 - my0) + 0.5)
            step_x = lx / ((lx + ly) / 2)
            sign = 1 if y1 > y0 else -1
            step_y = sign * ly / ((lx + ly) / 2)
            x_index = np.round(np.arange(mx0, mx1, step_x)).astype(int)
            y_index = np.round(np.arange(my0, my1, step_y)).astype(int)
            if len(x_index) != len(y_index):
                print('lengths differ')
            for x, y in zip(x_index, y_index):
                if y1 > y0:
                    mat_ff[[x], [y]] += pid
                    mat_ff[[y], [x]] += pid
                else:
                    mat_rc[[x], [y]] += pid
                    mat_rc[[y], [x]] += pid
    return [mat_ff, mat_rc]


def get_sequence_info(seq_file):
    """ Reads fasta sequences and
    return  number of sequences and total length"""
    with open(seq_file) as f:
        n_seq = 0
        l_seq = 0
        for line in f:
            if line.strip()[0] == ">":
                n_seq += 1
                continue
            else:
                l_seq += len(line.strip())
    return n_seq, l_seq


def blast(
        seq_file, word_lengths=(9, 13, 17, 23), ):
    """run series of blast with different word lengths """
    blast_params = "-outfmt 6 -ungapped -max_target_seqs 9999999 -dust no " \
                   "-perc_identity 95"

    subprocess.check_call(["makeblastdb", "-in", seq_file, "-dbtype", "nucl"])
    pp = list()
    out = list()

    for w in word_lengths:
        out.append(tempfile.NamedTemporaryFile(delete=False, mode='w').name)
        cmd = f'blastn -task megablast -db {seq_file} -query {seq_file} -out ' \
              f'{out[-1]} -word_size {w} {blast_params}'.split()

        pp.append(subprocess.Popen(cmd))
    for p in pp:
        p.wait()
    return out


class InvalidInputError(Exception):
    """ exception for wrong input file"""
    pass


def main():
    """
input - sequence in FASTA format
run blast against itself and returns blast tables converted to density matrix
    """
    description = "convert table with blast hits to matrix of densities"  # arguments
    # parsing
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-i', '--input', type=str, help="input fasta sequence", required=True
        )
    parser.add_argument(
        '-L', '--sequence_length', type=int, help="length if sequence ", required=True
        )
    parser.add_argument(
        '-W', '--array_width', type=int, required=True, help="required dimension of "
                                                             "output matrix"
        )
    parser.add_argument('-o', '--output', type=str, help="output file name prefix")
    args = parser.parse_args()
    number_of_seqs, length_of_seqs = get_sequence_info(args.input)
    print(f'input file with {number_of_seqs} sequence(s), total length {length_of_seqs}')
    if number_of_seqs != 1:
        raise InvalidInputError('input FASTA must contain single sequence')
    word_lengths = (9, 11)
    blast_tables = blast(args.input, word_lengths=word_lengths)
    for idx, val in enumerate(blast_tables):
        matrix_ff, matrix_rc = blast2density(
            w=args.array_width, length_of_seq=length_of_seqs, input_blast=val
            )
        # noinspection PyTypeChecker
        np.savetxt(
            fname=f'{args.output}_w{word_lengths[idx]}_F.csv', X=matrix_ff, delimiter="\t"
            )
        # noinspection PyTypeChecker
        np.savetxt(
            fname=f'{args.output}_w{word_lengths[idx]}_F.csv', X=matrix_rc, delimiter="\t"
            )


if __name__ == '__main__':
    main()
