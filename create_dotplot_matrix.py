#!/usr/bin/env python
""" convert alignment information from blast table to dotplot like density matrix"""
import itertools
import math
import multiprocessing
import os
import numpy as np
import csv
import argparse
import subprocess
import tempfile


def blast2density_and_export(w, length_of_seq, input_blast, word, pid, output_prefix):
    """ wrapper - run blast to density and write output"""
    matrix_ff, matrix_rc = blast2density(w, length_of_seq, input_blast)
    np.savetxt(
        fname=f'{output_prefix}_pid{pid}_w{word}_F.csv', X=matrix_ff, delimiter="\t"
        )
    np.savetxt(
        fname=f'{output_prefix}_pid{pid}_w{word}_R.csv', X=matrix_rc, delimiter="\t"
        )


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
                    mat_ff[[x], [y]] += 1
                    mat_ff[[y], [x]] += 1
                else:
                    mat_rc[[x], [y]] += 1
                    mat_rc[[y], [x]] += 1
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
        seq_file, word_lengths=(9, 13, 17, 23), perc_identity=(90, 95, 100)
        ):
    """run series of blast with different word lengths """
    blast_params = "-outfmt 6 -ungapped -max_target_seqs 9999999 -dust no "

    subprocess.check_call(["makeblastdb", "-in", seq_file, "-dbtype", "nucl"])
    pp = list()
    out = list()
    print('running blast')
    for w in word_lengths:
        for pid in perc_identity:
            out.append(tempfile.NamedTemporaryFile(delete=False, mode='w').name)
            cmd = f'blastn -task megablast -db {seq_file} -query {seq_file} -out ' \
                  f'{out[-1]} -word_size {w} {blast_params} -perc_identity ' \
                  f'{pid}'.split()
            print(cmd)
            pp.append(subprocess.Popen(cmd))

    for p in pp:
        p.wait()
    print('done')
    return out


class InvalidInputError(Exception):
    """ exception for wrong input file"""
    pass


def ccw(a, b, c):
    """ helper function fo intersection"""
    return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])


# Return true if line segments AB and CD intersect
def intersect(ax, ay, bx, by):
    """ return true is segments intersect"""
    return ccw(ax, bx, by) != ccw(ay, bx, by) and ccw(ax, ay, bx) != ccw(ax, ay, by)


def det(a, b):
    """ helper function fo intersection"""
    return a[0] * b[1] - a[1] * b[0]


def segment_intersection_relative(line1, line2):
    """ assumes that segments are always orthogonal"""
    if intersect(line1[0], line1[1], line2[0], line2[1]):
        d1 = math.dist(line1[0], line2[0])
        d2 = math.dist(line1[1], line2[0])
        return d1/(d1+d2)
    else:
        return False


def segment_intersection(line1, line2):
    """return position of intersection or false"""
    if intersect(line1[0], line1[1], line2[0], line2[1]):
        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
        div = det(xdiff, ydiff)
        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return x, y
    else:
        return False

def blast_orthogonal_intersects(input_blast, orth_segs):
    w = orth_segs[0][0][0]
    orth_hits = np.zeros([w + 1,len(orth_segs) + 1])
    print(w)
    print(len(orth_segs))
    with open(input_blast, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t', )
        for line in csvreader:
            if line[0][0] == "#":
                continue
            # coordinates of alignment
            y0 = float(line[8])
            y1 = float(line[9])
            # only plus-plu alignments
            if y1 < y0:
                continue
            x0 = float(line[6])
            x1 = float(line[7])
            for index, a in enumerate(orth_segs):
                #intersection = segment_intersection(a, ((x0,y0),(x1,y1)))
                intersection = segment_intersection(a, ((x0,y0),(x1,y1)))
                if intersection:
                    p = round(a[0][0] - intersection[0])
                    orth_hits[[p],[index]] = 1
    np.savetxt("tmp/test_orth3.csv", X=orth_hits)

def extract_diagonal_matrix_from_blast(input_blast, output, length_of_seq, w=5000,
                                       cfactor=100):
    """ return  matrix parallel with diagonal, max with of matrix is w"""

    diag_mat=np.full([w, length_of_seq], False)
    print('parsing blast')
    with open(input_blast) as f:
        csvreader = csv.reader(f, delimiter='\t', )
        for line in csvreader:
            if line[0][0] == "#":
                continue
            # coordinates of alignment
            y0 = int(line[8])
            y1 = int(line[9])
            # only plus-plu alignments
            if y1 < y0:
                continue
            x0 = int(line[6])
            x1 = int(line[7])
            d_dist = abs(x0 - y0)
            if d_dist >= 5000:
                continue
            if x0 > y0:
                d_position = range(x0, x1)
            else:
                d_position = range(y0, y1)
            diag_mat[[d_dist], d_position] = True
    # compress matrix:
    diag_mat_comp=np.full([w, int(length_of_seq/cfactor) + 1],0)
    print('compressing table')
    for i in range(int(length_of_seq/cfactor)):
        index=range(i * cfactor, (i + 1)* cfactor)
        diag_mat_comp[:,i] = np.sum(diag_mat[:,index], axis=1)
    np.savetxt(output, diag_mat_comp, fmt="%i")


def int_list(string):
    """
custom type of argparse
    :param string: int separated by comma (1,2,4,6)
    :return: list of integers
    """
    x = string.strip().split(",")
    xint = [int(i) for i in x]
    return xint


def main():
    """
input - sequence in FASTA format
run blast against itself and returns blast tables converted to density matrix
    """
    description = "Run blast on sequence against itself and blast hits to dotplot " \
                  "matrix of similarities"
    # arguments parsing
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-i', '--input', type=str, help="input fasta sequence", required=True
        )
    parser.add_argument(
        '-W', '--word_lengths', type=int_list, help="input word size for blast separated "
                                                    "by "
                                                    "comma", required=False, default="21,"
                                                                                     "31,41,51"
        )

    parser.add_argument(
        '-P', '--perc_identity', type=int_list, help="percentage identity thresholds"
                                                     " separated by comma",
                                                     required=False, default="90,95,100"
        )

    parser.add_argument(
        '-M', '--matrix_width', type=int, default=3000, required=False, help="required " \
                                                                        "dimension "
                                                                   "of "
                                                              "output matrix"
        )

    parser.add_argument(
        '-S', '--scale_coef', type=float, default=None, required=False,
        help="Matrix_width is length_of_seq/scale_coef"
        )

    parser.add_argument(
        '-o', '--output', type=str, help="output directory, will be "
                                         "created"
        )

    args = parser.parse_args()
    os.mkdir(args.output)
    output_prefix = f'{args.output}/density'
    param_file = f'{args.output}/params.txt'
    number_of_seqs, length_of_seqs = get_sequence_info(args.input)
    print(f'input file with {number_of_seqs} sequence(s), total length {length_of_seqs}')
    if number_of_seqs != 1:
        raise InvalidInputError('input FASTA must contain single sequence')
    word_lengths = args.word_lengths
    perc_identity = args.perc_identity

    blast_tables = blast(
        args.input, word_lengths=word_lengths, perc_identity=perc_identity
        )

    pid_word = itertools.product(word_lengths, perc_identity)

    params_list = []
    step_size = 100
    win_size = 15000

    diagonal_matrix_file = f'{args.output}/diag_matrix.csv'
    # run it only on one word_length/pid combination
    extract_diagonal_matrix_from_blast(blast_tables[0], length_of_seq=length_of_seqs,
                                       w=win_size, output = diagonal_matrix_file)

    if args.scale_coef:
        args.matrix_width = int(length_of_seqs / args.scale_coef)
    # calculate doplot like matrix
    for val in blast_tables:
        word, pid = next(pid_word)
        params_list.append(
            [args.matrix_width, length_of_seqs, val, word, pid, output_prefix]
            )

    with multiprocessing.Pool(20) as pool:
        pool.starmap(blast2density_and_export, params_list)
    [os.remove(i) for i in blast_tables]

    # calculate periodicity spectra - R script
    print('running get_periodicity_spectra')
    rscript = os.path.dirname(os.path.abspath(__file__)) + "/get_periodicity_spectra.R"
    # get periodicity spectra directly from sequence
    subprocess.check_call([rscript, '-i', args.input, '-t', 'FASTA', '-o', args.output +
                           "/periodicity_fft.RDS"])
    # get periodicity spectra from blast output
    subprocess.check_call(
        [rscript, '-i', diagonal_matrix_file, '-t', 'MATRIX','-o',
         args.output + "/periodicity_blast_fft.RDS"]
        )
    # export parameters:
    # note this param file is used in follow up script
    # check for compatibility before changing it
    with open(param_file, mode='w') as f:
        print(f'sequence_length : {length_of_seqs}', file=f)
        print(f'word_lengths    : {args.word_lengths}', file=f)
        print(f'matrix width    : {args.matrix_width}', file=f)
        print(f'input file      : {args.input}', file=f)


if __name__ == '__main__':
    main()
