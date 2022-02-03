
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2022/02/01
# Version     : 1.0
# © Copyright : 2022 Dominique Zeise
# =============================================================================
"""
Multiple Sequence Alignment (MSA)
Global Alignment Algorithm


This module uses dynamic programming to find the optimal
alignment of multiple sequences/ strings of similar length.
Its the optimized version of the algorithm.
"""
# =============================================================================


import numpy as np

MATCH = 1
MISMATCH = -1
GAP = -2


def fill_matches(score_matrix, ls1, ls2, ls3, seq1, seq2, seq3):
    """creates a substitution matrix inside
    initialized matrix.
    Starts with second row and column.

    Args:
        score_matrix (np matrix): initialized matrix
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt
        ls3 (int): len third sequence to alignt
        seq1 (str): first sequence to alignt
        seq2 (str): second sequence to alignt
        seq2 (str): third sequence to alignt

    Returns:
        numpy matrix: substitution matrix inside
    initialized matrix
    """
    for z in range(ls3):
        for y in range(ls2):
            for x in range(ls1):
                if ((seq1[x] == seq2[y]) or (seq1[x] == seq3[z]) or (seq2[y] == seq3[z])):
                    score_matrix[z][y][x] = MATCH
                else:
                    score_matrix[z][y][x] = MISMATCH
    return score_matrix


def initialize_matrix(score_matrix, ls1, ls2, ls3):
    """fillst first row and column with
    given penalties

    Args:
        score_matrix (numpy matrix): 0 filled
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt
        ls3 (int): len third sequence to alignt

    Returns:
        numpy matrix: np matrix with initialized
        first row and column
    """
    for x in range(ls1):
        score_matrix[0][0][x] = x * 2*GAP
    for y in range(ls2):
        score_matrix[0][y][0] = y * 2*GAP
    for z in range(ls3):
        score_matrix[z][0][0] = z * 2*GAP
    return score_matrix


def recalculate_scorematrix(score_matrix, ls1, ls2, ls3):
    """creates a scorematrix with help of the
    substitution matrix and rewards/ penalties

    Args:
        score_matrix (np matrix): substitution matrix inside
        initialized matrix
        ls1 (int): len first sequence to alignt
        ls2 (int): len second sequence to alignt
        ls3 (int): len third sequence to alignt

    Returns:
        numpy matrix: scorematrix
    """
    # start at second row and column -> skips initialized cells
    for z in range(ls3):
        for y in range(ls2):
            for x in range(ls1):
                coords = (x,y,z)
                # surfaces
                if coords.count(0) == 1:
                    # front surface
                    if z == 0:
                        score_matrix[z][y][x] = max(
                            score_matrix[z][y-1][x-1] + GAP,
                            score_matrix[z][y-1][x] + 2*GAP,
                            score_matrix[z][y][x-1] + 2*GAP
                        )
                    # top surface
                    elif y == 0:
                        score_matrix[z][y][x] = max(
                            score_matrix[z-1][y][x-1] + GAP,
                            score_matrix[z-1][y][x] + 2*GAP,
                            score_matrix[z][y][x-1] + 2*GAP
                        )
                    # left side surface
                    else:
                        score_matrix[z][y][x] = max(
                            score_matrix[z-1][y-1][x] + GAP,
                            score_matrix[z-1][y][x] + 2*GAP,
                            score_matrix[z][y-1][x] + 2*GAP
                        )

                # border behaviour
                elif coords.count(0) == 2:
                    # front top border -> go to left
                    if z == y == 0:
                        score_matrix[z][y][x] = score_matrix[z][y][x-1] + 2*GAP
                    # front left border -> go up
                    elif z == x == 0:
                        score_matrix[z][y][x] = score_matrix[z][y-1][x] + 2*GAP
                    # top left border -> go to front
                    elif z == x == 0:
                        score_matrix[z][y][x] = score_matrix[z-1][y][x] + 2*GAP

                # start corrdinates
                elif coords.count(0) == 3:
                    continue

                # inside cube
                else:
                    score_matrix[z][y][x] = max(
                        # diagonal (3d movement)
                        score_matrix[z-1][y-1][x-1] + score_matrix[z][y][x],
                        # diagonal neighbors (2d movement)
                        score_matrix[z-1][y][x-1] + GAP,
                        score_matrix[z-1][y-1][x] + GAP,
                        score_matrix[z][y-1][x-1] + GAP,
                        # left faces (1d movement)
                        score_matrix[z][y][x-1] + 2*GAP,
                        score_matrix[z-1][y][x] + 2*GAP,
                        score_matrix[z][y-1][x] + 2*GAP,
                    )

    return score_matrix


def traceback(score_matrix, ls1, ls2, ls3, seq1, seq2, seq3):
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
    seq1_new, seq2_new, seq3_new = [], [], []
    x = ls1-1
    y = ls2-1
    z = ls3-1

    while x > 0 and y > 0 and z > 0:
        # diagonal (3d movement)
        xyz_dia = score_matrix[z-1][y-1][x-1]
        # diagonal neighbors (2d movement)
        xy_dia = score_matrix[z][y-1][x-1]
        xz_dia = score_matrix[z-1][y][x-1]
        yz_dia = score_matrix[z-1][y-1][x]
        # left faces (1d movement)
        left = score_matrix[z][y][x-1]
        up = score_matrix[z][y-1][x]
        front = score_matrix[z-1][y][x]

        max_points = max(xyz_dia, xz_dia, yz_dia, xy_dia, left, front, up)

        if xyz_dia == max_points:
            seq1_new.append(seq1[x-1])
            seq2_new.append(seq2[y-1])
            seq3_new.append(seq3[z-1])
            x -= 1
            y -= 1
            z -= 1
        elif yz_dia == max_points:
            seq1_new.append("_")
            seq2_new.append(seq2[y-1])
            seq3_new.append(seq3[z-1])
            y -= 1
            z -= 1
        elif xy_dia == max_points:
            seq1_new.append(seq1[x-1])
            seq2_new.append(seq2[y-1])
            seq3_new.append("_")
            x -= 1
            y -= 1
        elif xz_dia == max_points:
            seq1_new.append(seq1[x-1])
            seq2_new.append("_")
            seq3_new.append(seq3[z-1])
            x -= 1
            z -= 1

        elif front == max_points:
            seq1_new.append("_")
            seq2_new.append("_")
            seq3_new.append(seq3[z-1])
            z -= 1

        elif up == max_points:
            seq1_new.append("_")
            seq2_new.append(seq2[y-1])
            seq3_new.append("_")
            y -= 1
        # left
        else:
            seq1_new.append(seq1[x-1])
            seq2_new.append("_")
            seq3_new.append("_")
            x -= 1
    return seq1_new, seq2_new, seq3_new


def output(seq1_new, seq2_new, seq3_new):
    """generates and prints an alignment output string
    which shows matches (|) and mismatches (*)

    Args:
        seq1_new (str): redesigned first sequence
        seq2_new (str): redesigned second sequence

    Returns:
        str: alignment output string which
        consists of "|" and "*" symbols
    """
    alignment1_2 = []
    alignment1_3 = []
    alignment2_3 = []

    for n1, n2, n3 in zip(seq1_new, seq2_new, seq3_new):
        alignment1_2.append("|") if (n1 == n2) and (
            n1 != "_") else alignment1_2.append("*")
        alignment1_3.append("|") if (n1 == n3) and (
            n1 != "_") else alignment1_3.append("*")
        alignment2_3.append("|") if (n2 == n3) and (
            n2 != "_") else alignment2_3.append("*")
    seq1_output = ''.join(seq1_new[:-1])[::-1]
    seq2_output = ''.join(seq2_new[:-1])[::-1]
    seq3_output = ''.join(seq3_new[:-1])[::-1]
    alignment_output1_2 = ''.join(alignment1_2[:-1])[::-1]
    alignment_output1_3 = ''.join(alignment1_3[:-1])[::-1]
    alignment_output2_3 = ''.join(alignment2_3[:-1])[::-1]

    print("Sequence 1 vs Sequence 2")
    print(seq1_output)
    print(alignment_output1_2)
    print(seq2_output, "\n")

    print("Sequence 1 vs Sequence 3")
    print(seq1_output)
    print(alignment_output1_3)
    print(seq3_output, "\n")

    print("Sequence 2 vs Sequence 3")
    print(seq2_output)
    print(alignment_output2_3)
    print(seq3_output)

    return alignment_output1_2, alignment_output1_3, alignment_output2_3


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


def main():
    # change sequences below for your needs!
    seq1 = ",AA,"
    seq2 = ",AT,"
    seq3 = ",AA,"

    ls1 = len(seq1)
    ls2 = len(seq2)
    ls3 = len(seq3)
    score_matrix = np.zeros((ls3, ls2, ls1), dtype="int16")
    score_matrix = fill_matches(score_matrix, ls1, ls2, ls3, seq1, seq2, seq3)
    score_matrix = initialize_matrix(score_matrix, ls1, ls2, ls3)
    score_matrix = recalculate_scorematrix(score_matrix, ls1, ls2, ls3)
    seq1_new, seq2_new, seq3_new = traceback(score_matrix, ls1, ls2, ls3, seq1, seq2, seq3)
    alignment1_2, alignment1_3, alignment2_3 = output(
        seq1_new, seq2_new, seq3_new)
    for alignment in (alignment1_2, alignment1_3, alignment2_3):
        sim_tup = calc_similarity(alignment)
        print(f"\nSimilarity: {round(sim_tup[1]*100, 2)}%")


if __name__ == '__main__':
    main()
