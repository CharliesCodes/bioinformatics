
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2021/11/03
# Version     : 1.0
# © Copyright : 2021 Dominique Zeise
# =============================================================================
"""
Free-Shift Alignment/ Semiglobales Alignment/ End-Gap Free Alignment


This module uses dynamic programming to find the optimal
alignment of two sequences/ strings of different length.
It ignores overlaps (prefixes/ suffixes).
"""
# =============================================================================

import numpy as np

MATCH = 1
MISMATCH = -10
GAP = -10


def fill_matches(score_matrix, ls1, ls2, seq1, seq2):
    """creates a substitution matrix.
    Starts with second row and column.

    Args:
        score_matrix (np matrix): initialized matrix
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt
        seq1 (str): first sequence to alignt
        seq2 (str): second sequence to alignt

    Returns:
        numpy matrix: substitution matrix inside
    initialized matrix
    """
    for y in range(1, ls2):
        for x in range(1, ls1):
            score_matrix[y][x] = MATCH if seq1[x] == seq2[y] else MISMATCH
    return score_matrix


def recalculate_scorematrix(score_matrix, ls1, ls2):
    """creates a scorematrix with help of the
    substitution matrix and rewards/ penalties

    Args:
        score_matrix (np matrix): substitution matrix inside
        initialized matrix
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt

    Returns:
        numpy matrix: scorematrix
    """
    # start at second row and column -> skips initialized cells
    for y in range(1, ls2):
        for x in range(1, ls1):
            score_matrix[y][x] = max(
                score_matrix[y-1][x-1] + score_matrix[y][x],
                score_matrix[y-1][x] + GAP,
                score_matrix[y][x-1] + GAP
            )
    return score_matrix


def get_max_from_border(score_matrix):
    """find coordinates of highest score in the
    first row and column of the matrix

    Args:
        score_matrix (numpy matrix): scorematrix

    Returns:
        tuple: (y, x) Coordinates of max score
        returns the first occuring if score occures multiple times
    """
    matrix_shape = score_matrix.shape
    y = matrix_shape[0]-1
    x = matrix_shape[1]-1

    last_column, last_row = score_matrix[:, -1], score_matrix[-1]
    last_col_max, last_row_max = np.max(last_column), np.max(last_row)

    if last_col_max >= last_row_max:
        result = (np.argmax(last_column), x)
    else:
        result = (y, np.argmax(last_row))
    return result


def traceback(score_matrix, seq1, seq2, max_border_coords):
    """traceback from coordinates of highest score in
    last row or last column to upper left corner.
    Finds best scoring path via greedy-algorithm.

    Args:
        score_matrix numpy matrix: scorematrix
        seq1 (str): first sequence to alignt
        seq2 (str): second sequence to alignt
        max_border_coords (tuple): (y, x) Coordinates of max score

    Returns:
        seq1_new (str): redesigned first sequence
        seq2_new (str): redesigned second sequence
    """
    seq1_new, seq2_new = [], []
    y, x = max_border_coords[0], max_border_coords[1]
    seq1_new.append(seq1[x])
    seq2_new.append(seq2[y])

    while x > 0 and y > 0:
        dia = score_matrix[y-1][x-1]
        up = score_matrix[y-1][x]
        left = score_matrix[y][x-1]
        max_points = max(dia, up, left)

        # dia
        if dia == max_points:
            seq1_new.append(seq1[x-1])
            seq2_new.append(seq2[y-1])
            x -= 1
            y -= 1
        # up
        elif up == max_points:
            seq1_new.append("_")
            seq2_new.append(seq2[y-1])
            y += -1
        # left
        else:
            seq2_new.append("_")
            seq1_new.append(seq1[x-1])
            x += -1

    # left border
    while y > 0:
        seq1_new.append("_")
        seq2_new.append(seq2[y-1])
        y += -1
    # top border
    while x > 0:
        seq2_new.append("_")
        seq1_new.append(seq1[x-1])
        x += -1
    return seq1_new, seq2_new


def add_overlap(ls1, ls2, seq1, seq2, seq1_new, seq2_new, max_border_coords):
    """merges original and new sequences together and
    appends underscore symbols at the shorter sequence as suffix

    Args:
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt
        seq1 (str): first sequence to alignt
        seq2 (str): second sequence to alignt
        seq1_new (str): redesigned first sequence
        seq2_new (str): redesigned second sequence
        max_border_coords (tuple): (y, x) Coordinates of max score

    Returns:
        str: redesigned sequences with suffixes
    """
    y, x = max_border_coords
    seq1_new = ''.join(seq1_new)[::-1]
    seq2_new = ''.join(seq2_new)[::-1]
    # add suffix overlap
    seq1_new = seq1_new + seq1[x+1:] + "_"*abs(ls2-1 - y)
    seq2_new = seq2_new + seq2[y+1:] + "_"*abs(ls1-1 - x)
    return seq1_new, seq2_new


def merge_sequences(seq1_new, seq2_new):
    """used for assembly only!
    merges both sequences to one at overlap

    Args:
        seq1_new (str): main sequence
        seq2_new (str): seq to append to

    Returns:
        seq1_new (str): merged main sequence
    """
    seq1_new = list(seq1_new)
    for x in range(len(seq1_new)):
        if seq1_new[x] == "_":
            seq1_new[x] = seq2_new[x]
        if seq1_new[x] == ",":
            del seq1_new[x]
    seq1_new = ''.join(seq1_new)
    return seq1_new


def output(seq1_new, seq2_new):
    """generates and prints an alignment output string
    which shows matches (|) and mismatches (*)

    Args:
        seq1_new (str): redesigned first sequence
        seq2_new (str): redesigned second sequence

    Returns:
        str: alignment output string which
        consists of "|" and "*" symbols
    """
    alignment = []
    for n1, n2 in zip(seq1_new, seq2_new):
        alignment.append("|") if n1 == n2 else alignment.append("*")
    alignment_output = ''.join(alignment)
    return alignment_output


def calc_similarity(alignment_output):
    """calculates absolute & relative match frequency
    of the alignment by counting the Pipe symbols
    in the alignment_output string

    Args:
        alignment_output (str): consists of "|" and "*" symbols
        for matches and mismatches in the alignment

    Returns:
        tuple: absolute & relative match frequency
    """
    abs_count = alignment_output.count("|")
    rel_count = abs_count/len(alignment_output)
    return (abs_count, rel_count)


def main(seq1='', seq2='', assembly=False):
    # change sequences below for your needs!
    if not (seq1 or seq2):
        seq1 = ",ATTAC"
        seq2 = ",ATT"
    ls1 = len(seq1)
    ls2 = len(seq2)
    score_matrix = np.zeros(ls2 * ls1, dtype="int16").reshape((ls2, ls1))
    score_matrix = fill_matches(score_matrix, ls1, ls2, seq1, seq2)
    score_matrix = recalculate_scorematrix(score_matrix, ls1, ls2)
    max_border_coords = get_max_from_border(score_matrix)
    seq1_new, seq2_new = traceback(
        score_matrix, seq1, seq2, max_border_coords)
    seq1_new, seq2_new = add_overlap(
        ls1, ls2, seq1, seq2, seq1_new, seq2_new, max_border_coords)
    # remove Commas
    seq1_new = seq1_new.replace(',', '')
    seq2_new = seq2_new.replace(',', '')

    alignment_output = output(seq1_new, seq2_new)
    sim_tup = calc_similarity(alignment_output)

    if assembly:
        seq1_new = merge_sequences(seq1_new, seq2_new)
        return sim_tup, seq1_new
    else:
        print(seq1_new)
        print(alignment_output)
        print(seq2_new)
        print(f"\nAehnlichkeit: {round(sim_tup[1]*100, 2)}%")


if __name__ == '__main__':
    main()
