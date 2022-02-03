
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2021/11/03
# Version     : 2.0
# © Copyright : 2021 Dominique Zeise
# =============================================================================
"""
Needleman-Wunsch Algorithm
Global Alignment Algorithm


This module uses dynamic programming to find the optimal
alignment of two sequences/ strings of similar length.
Its the optimized version of the algorithm.
"""
# =============================================================================

import numpy as np

MATCH = 1
MISMATCH = -1
GAP = -2


def initialize_matrix(score_matrix, ls1, ls2):
    """fillst first row and column with
    given penalties

    Args:
        score_matrix (numpy matrix): 0 filled
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt

    Returns:
        numpy matrix: np matrix with initialized
        first row and column
    """
    for x in range(ls1):
        score_matrix[0][x] = x * GAP
    for y in range(ls2):
        score_matrix[y][0] = y * GAP
    return score_matrix


def fill_matches(score_matrix, ls1, ls2, seq1, seq2):
    """creates a substitution matrix inside
    initialized matrix.
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


def traceback(score_matrix, ls1, ls2, seq1, seq2):
    """traceback from bottom right corner to
    upper left corner. Finds best scoring path
    via greedy-algorithm.

    Args:
        score_matrix numpy matrix: scorematrix
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt
        seq1 (str): first sequence to alignt
        seq2 (str): second sequence to alignt

    Returns:
        seq1_new (str): redesigned first sequence
        seq2_new (str): redesigned second sequence
    """
    seq1_new, seq2_new = [], []
    x = ls1-1
    y = ls2-1
    while x>0 and y>0:
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
            y -= 1
        # left
        else:
            seq2_new.append("_")
            seq1_new.append(seq1[x-1])
            x -= 1
    # left border
    while y > 0:
        seq1_new.append("_")
        seq2_new.append(seq2[y-1])
        y -= 1
    # top border
    while x > 0:
        seq2_new.append("_")
        seq1_new.append(seq1[x-1])
        x -= 1
    return seq1_new, seq2_new


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
    seq1_output = ''.join(seq1_new[:-1])[::-1]
    seq2_output = ''.join(seq2_new[:-1])[::-1]
    alignment_output = ''.join(alignment[:-1])[::-1]

    print(seq1_output)
    print(alignment_output)
    print(seq2_output)
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


def main(seq1='', seq2=''):
    # change sequences below for your needs!
    if not (seq1 or seq2):
        seq1 = ",AATGC,"
        seq2 = ",ATGC,"
    ls1 = len(seq1)
    ls2 = len(seq2)
    score_matrix = np.zeros((ls2, ls1), dtype="int16")
    score_matrix = initialize_matrix(score_matrix, ls1, ls2)
    score_matrix = fill_matches(score_matrix, ls1, ls2, seq1, seq2)
    score_matrix = recalculate_scorematrix(score_matrix, ls1, ls2)
    seq1_new, seq2_new = traceback(score_matrix, ls1, ls2, seq1, seq2)
    alignment_output = output(seq1_new, seq2_new)
    sim_tup = calc_similarity(alignment_output)
    print(f"\nSimilarity: {round(sim_tup[1]*100, 2)}%")


if __name__ == '__main__':
    main()
