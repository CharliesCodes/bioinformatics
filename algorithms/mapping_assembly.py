
# =============================================================================
# Created By  : Dominique Zeise
# GitHub      : https://github.com/CharliesCodes
# Created Date: 2021/12/22
# Version     : 1.0
# © Copyright : 2021 Dominique Zeise
# =============================================================================
'''Mapping Assembly


This module uses dynamic programming to restore an origin-sequence from
sequenceparts.
Dependencies: Local Alignment Module'''
# =============================================================================
try:
    import free_shift_alignment as fsa
except ImportError:
    print("Free-Shift Alignment not found!\nPlease copy file 'free_shift_alignment.py' into the same directory!\nVisit my GitHub page to download it.")

import random


def generate_origin(origin_len=100):
    origin = ''.join(random.choice(["A", "T", "C", "G"]) for x in range(origin_len))
    return origin, origin_len


def cut_origin_to_seqparts(origin, minl, maxl, sets=20):
    """This func generates random sets of sequenceparts from origin
    and wraps them into one big, sorted list

    Args:
        origin (str): original sequence to cut
        minl (int): min len for sequenceparts
        maxl (int): max len for sequenceparts
        sets (int, optional): number of sets from original. Defaults to 10.

    Returns:
        [list]: len decreasing sorted  list of substrings with random parts from origin
    """
    seqparts = []
    for _ in range(sets):
        x = 0
        while x < len(origin):
            rand = random.randint(minl, maxl)
            if (len(origin)-x) <= minl:
                seqparts.append(origin[x::])
            else:
                seqparts.append(origin[x:x+rand])
            x += rand
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


def mapping(seqparts, minl, origin_len):
    seqparts = ["," + seq for seq in seqparts]
    # set longest (first) list element as main sequence and remove it from list
    main_seq = seqparts.pop(0)
    while seqparts and len(main_seq) <= origin_len:
        fragments = []
        for seq2 in seqparts:
            # use imported smith-waterman algorithm without output function
            sim_tup, main_seq_new = fsa.main(main_seq, seq2, True)
            if 0.9*len(seq2) >= sim_tup[0] > 0.2*len(seq2):
                fragments.append((sim_tup[0], main_seq_new))
        if fragments:
            sorted_frags = sorted(fragments, key=lambda tup: tup[0], reverse=True)
            main_seq = sorted_frags[0][1]
            main_seq = "," + main_seq
        seqparts.pop(0)
    return main_seq


def main():
    origin, origin_len = generate_origin()
    minl = round(0.10*origin_len)
    maxl = round(0.20*origin_len)
    seqparts = cut_origin_to_seqparts(origin, minl, maxl, sets=20)
    seqparts = get_unique_seqparts(seqparts)
    assembly_sequence = mapping(seqparts, minl, origin_len)
    # global alignment here:
    origin = "," + origin +","
    assembly_sequence = assembly_sequence + ","
    print("Ursprungssequenz")
    print("Assemblierte Sequenz\n")
    fsa.main(origin, assembly_sequence)


if __name__ == '__main__':
    main()
