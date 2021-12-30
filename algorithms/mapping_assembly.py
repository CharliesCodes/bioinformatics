
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
    import smith_waterman as sw
except ImportError:
    print("Smith-Waterman Algorithm not found!\nPlease copy file 'smith_waterman.py' into the same directory!\nVisit my GitHub page to download it.")

import random


def generate_origin(origin_len=100):
    origin = ''.join(random.choice(["A", "T", "C", "G"]) for x in range(origin_len))
    return origin, origin_len


def cut_origin_to_seqparts(origin, minl, maxl, sets=10):
    """This func generates random sets of sequenceparts from origin
    and wraps them into one big, sorted list

    Args:
        origin (str): original sequence to cut
        minl (int): min len for sequenceparts
        maxl (int): max len for sequenceparts
        sets (int, optional): number of sets from original. Defaults to 10.

    Returns:
        [list]: list of substrings with random parts from origin
    """
    seqparts = []
    for _ in range(sets):
        x = 0
        while x < len(origin):
            rand = random.randint(minl, maxl)
            seqparts.append(origin[x:x+rand])
            if (len(origin)-x) <= minl:
                seqparts.append(origin[x::])
            x += rand
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
    for substring in seqparts:
        if not any(substring in string for string in unique_seqparts):
            unique_seqparts.append(substring)
    return unique_seqparts


def test(seqparts):
    #! SOME CODESMELL is here
    seqparts = ["," + seq + "," for seq in seqparts]
    main_seq = seqparts[0]
    seqparts.pop(0)
    while seqparts:
        fragments = []
        for seq2 in seqparts:
            sim_tup, seq1_new = sw.main(main_seq, seq2)
            fragments.append((sim_tup[0], seq1_new))
        sorted_frags = sorted(fragments, key=lambda tup: tup[0], reverse=True)
        main_seq = sorted_frags[0][1]
        print(sorted_frags)


def main():
    origin, origin_len = generate_origin()
    minl = round(0.01*origin_len)
    maxl = round(0.10*origin_len)
    seqparts = cut_origin_to_seqparts(origin, minl, maxl, sets=5)
    seqparts = get_unique_seqparts(seqparts)
    test(seqparts)





if __name__ == '__main__':
    main()
