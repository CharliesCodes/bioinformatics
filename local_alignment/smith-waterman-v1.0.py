
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2021/11/03
# Version     : 1.0
# © Copyright : 2021 Dominique Zeise
# =============================================================================
"""
Smith-Waterman Algorithm


This module uses dynamic programming to finding the optimal
alignment of two sequences/ strings of different lengths.
Its the optimization version of the algorithm.
"""
# =============================================================================

import numpy as np

MATCH = 1
MISMATCH = -1
GAP = -2


def fill_matches(score_matrix, seq1, seq2):
    for y in range(1, len(seq2)):
        for x in range(1, len(seq1)):
            score_matrix[y][x] = MATCH if seq1[x] == seq2[y] else MISMATCH
    return score_matrix


def recalculate_scorematrix(score_matrix, seq1, seq2):
    for y in range(1, len(seq2)):
        for x in range(1, len(seq1)):
            score_matrix[y][x] = max(
                0,
                score_matrix[y-1][x-1] + score_matrix[y][x],
                score_matrix[y-1][x] + GAP,
                score_matrix[y][x-1] + GAP
                )
    return score_matrix


def find_max_coordinates(score_matrix):
    result = np.where(score_matrix == np.amax(score_matrix))
    max_coords = (result[0][0].item(), result[1][0].item())
    return max_coords


def traceback(score_matrix, seq1, seq2, max_coords):
    seq1_new, seq2_new = [], []
    (y, x) = max_coords
    seq1_new.append(seq1[x])
    seq2_new.append(seq2[y])

    while score_matrix[y][x] != 0:
        dia = score_matrix[y-1][x-1]
        up = score_matrix[y-1][x]
        left = score_matrix[y][x-1]
        max_points = max(dia, up, left)
        score_matrix[y][x] = 0
        # dia
        if (x > 0 and y > 0 and dia == max_points):
            seq1_new.append(seq1[x-1])
            seq2_new.append(seq2[y-1])
            x -= 1
            y -= 1
        # up
        elif(y > 0 and up == max_points):
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
        alignment.append("|") if n1 == n2 else alignment.append("*")
    seq1_output = ''.join(seq1_new[:-1])[::-1]
    seq2_output = ''.join(seq2_new[:-1])[::-1]
    alignment_output = ''.join(alignment[:-1])[::-1]

    print(seq1_output)
    print(alignment_output)
    print(seq2_output)


def main():
    seq1 = ",ACGACGACTTACTT"
    seq2 = ",CGTCT"

    score_matrix = np.zeros((len(seq2)) * (len(seq1)),
                            dtype="int16").reshape((len(seq2), len(seq1)))

    score_matrix = fill_matches(score_matrix, seq1, seq2)
    score_matrix = recalculate_scorematrix(score_matrix, seq1, seq2)
    max_coords = find_max_coordinates(score_matrix)
    seq1_new, seq2_new = traceback(score_matrix, seq1, seq2, max_coords)
    del score_matrix
    output(seq1_new, seq2_new)


if __name__ == '__main__':
    main()
