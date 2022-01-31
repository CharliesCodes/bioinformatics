
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2021/12/22
# Version     : 1.0
# © Copyright : 2021 Dominique Zeise
# =============================================================================
'''De novo Assembly


This module uses dynamic programming to create a
de novo transcriptome without a reference genome.

Dependencies: Free-Shift Alignment Module'''
# =============================================================================

import random

try:
    import free_shift_alignment as fsa
except ImportError:
    print("Free-Shift Alignment not found!\nPlease copy file 'free_shift_alignment.py' into the same directory!\nVisit my GitHub page to download it.")


MINLEN = 5
MAXLEN = 15
SUBSEQ_NUM = 200
ASSEMBLY_MAXLEN = 100


def create_seqparts():
    """This func generates random sets of sequenceparts
    and wraps them into one big, sorted list

    Returns:
        [list]: len decreasing sorted list of substrings with random parts from origin
    """
    seqparts = []
    for _ in range(SUBSEQ_NUM):
        rand = random.randint(MINLEN, MAXLEN)
        subseq = ''.join(random.choice(
            ["A", "T", "C", "G"]) for x in range(rand))
        seqparts.append(subseq)
    # sorting not needed - just for the sake of clarity
    seqparts.sort(key=len)
    return seqparts


def get_unique_seqparts(seqparts):
    """This func generates a list of unique sequences

    Args:
        seqparts (list): strings of unordered, ununique sequence parts

    Returns:
        list: unique sequence parts - order: size decreasing
    """
    seqparts = list(set(seqparts))

    seqparts.sort(key=len, reverse=True)
    unique_seqparts = []
    # remove all sequences that are part of any other sequence
    for substring in seqparts:
        if not any(substring in string for string in unique_seqparts):
            unique_seqparts.append(substring)
    return unique_seqparts


def mapping(seqparts):
    """assembles the best fitting sequence into
    one big sequence. Requires Free-Shift-Alignment
    to find best fit. Starting with biggest Sequencepart,
    it itterates over all other parts (size decreasing order).
    Stores scores for fit in tuple and merges highest fitting
    sequence. Repeating the procedure untill seqparts list is
    empty or the ASSEMBLY_MAXLEN limit is reached

    Args:
        seqparts (list): strings with random generated sequences
        Order: decreasing length

    Returns:
        str: Assembled Sequence
    """
    seqparts = ["," + seq for seq in seqparts]
    # set longest (first) list element as main sequence and remove it from list
    main_seq = seqparts.pop(0)
    # assembly start sequence
    print("Startsequence:", main_seq)
    while seqparts and len(main_seq) <= ASSEMBLY_MAXLEN:
        fragments = []
        for parts_index, seq2 in enumerate(seqparts):
            # use imported smith-waterman algorithm without output function
            sim_tup, main_seq_new = fsa.main(main_seq, seq2, True)
            if 0.9*len(seq2) >= sim_tup[0] > 0.2*len(seq2):
                fragments.append((sim_tup[0], main_seq_new, parts_index))
        if fragments:
            sorted_frags = sorted(
                fragments, key=lambda tup: tup[0], reverse=True)
            main_seq = sorted_frags[0][1]
            main_seq = "," + main_seq
        try:
            del seqparts[sorted_frags[0][2]]
        except:
            print("An exception occurred! No fitting part found")

    return main_seq


def main():
    seqparts = create_seqparts()
    seqparts = get_unique_seqparts(seqparts)
    assembly_sequence = mapping(seqparts)
    print("\nAssemblierte Sequenz")
    print(assembly_sequence)


if __name__ == '__main__':
    main()
