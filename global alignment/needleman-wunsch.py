
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2021/11/03
# Version     : 1.0
# © Copyright : 2021 Dominique Zeise
# =============================================================================
"""
Needleman-Wunsch Algorithm


This module uses dynamic programming to finding the optimal
alignment of two sequences/ strings.
"""
# =============================================================================

import numpy as np

MATCH = 1
MISMATCH = -1
GAP = -2


def fill_match_matrix(match_matrix, seq1, seq2):
    for y in range(len(seq2)):
        for x in range(len(seq1)):
            if seq1[x] == seq2[y]:
                match_matrix[y][x] = MATCH
            else:
                match_matrix[y][x] = MISMATCH
    return match_matrix


def initialize_matrix(score_matrix, seq1, seq2):
    for x in range(len(seq1)):
        score_matrix[0][x] = x * GAP
    for y in range(len(seq2)):
        score_matrix[y][0] = y * GAP
    return score_matrix


def fill_scorematrix(score_matrix, match_matrix, seq1, seq2):
    for x in range(1, len(seq1)):
        for y in range(1, len(seq2)):
            score_matrix[y][x] = max(
                score_matrix[y-1][x-1] + match_matrix[y-1][x-1],
                score_matrix[y-1][x] + GAP,
                score_matrix[y][x-1] + GAP
                )
    return score_matrix


def traceback(score_matrix, match_matrix, seq1, seq2):
    seq1_new, seq2_new = [], []
    x = len(seq1)-1
    y = len(seq2)-1

    while x>0 and y>0:
        # dia
        if (x>0 and y>0 and score_matrix[y][x] == (score_matrix[y-1][x-1] +
                                                   match_matrix[y-1][x-1])):
            seq1_new.append(seq1[x-1])
            seq2_new.append(seq2[y-1])
            x -= 1
            y -= 1
        # up
        elif(y>0 and score_matrix[y][x] == score_matrix[y-1][x] + GAP):
            seq1_new.append("_")
            seq2_new.append(seq2[y-1])
            y += -1
        # left
        else:
            seq2_new.append("_")
            seq1_new.append(seq1[x-1])
            x += -1
    return seq1_new, seq2_new


def output(seq1_new, seq2_new):
    alignment = []

    for n1, n2 in zip(seq1_new, seq2_new):
        if n1 == n2: alignment.append("|")
        else: alignment.append("*")
    seq1_output = ''.join(seq1_new[:-1])[::-1]
    seq2_output = ''.join(seq2_new[:-1])[::-1]
    alignment_output = ''.join(alignment[:-1])[::-1]

    print(seq1_output)
    print(alignment_output)
    print(seq2_output)


def main():

    seq1 = ",ATGCT,"
    seq2 = ",AGCT,"

    score_matrix = np.zeros((len(seq2)) * (len(seq1)),
                            dtype="int16").reshape((len(seq2), len(seq1)))
    match_matrix = np.zeros((len(seq2)) * (len(seq1)),
                            dtype="int8").reshape((len(seq2), len(seq1)))

    match_matrix = fill_match_matrix(match_matrix, seq1, seq2)
    score_matrix = initialize_matrix(score_matrix, seq1, seq2)
    score_matrix = fill_scorematrix(score_matrix, match_matrix, seq1, seq2)
    seq1_new, seq2_new = traceback(score_matrix, match_matrix, seq1, seq2)
    del match_matrix, score_matrix
    output(seq1_new, seq2_new)


if __name__ == '__main__':
    main()
