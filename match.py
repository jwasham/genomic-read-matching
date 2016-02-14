"""
Tool for finding subsequences (called reads) in a larger DNA sequence
even though subsequence may contain a single nucleotide error
"""

def get_sequence():
    """
    Returns this sequence:
    Homo sapiens ABO gene, partial cds, promoter region, allele: B3
    http://www.ncbi.nlm.nih.gov/nuccore/LC068776.1
    :return:
    """

    references = [
        'ggccgcctcc cgcgcccctc tgtcccctcc cgtgttcggc ctcgggaagt cggggcggcg',
        'ggcggcgcgg gccgggaggg gtcgcctcgg gctcaccccg ccccagggcc gccgggcgga',
        'aggcggaggc cgagaccaga cgcggagcca tggccgaggt gttgcggacg ctggccg',
    ]

    ref_string = ''.join(references)

    return ref_string.replace(' ', '')


def build_sequence_index(sequence, key_length):
    """
    Builds a dictionary of (key_length)-char sequence -> position in larger sequence
    :param sequence:
    :param key_length:
    :return:
    """

    index = {}

    for c in range(len(sequence) - key_length + 1):
        seq_fragment = sequence[c:c + key_length]
        index[seq_fragment] = c + 1

    return index


def get_position(part, key_length, index):
    """
    Find the position of given short sequence in the index

    Since we can expect at most one error in the subsequence, if
    we don't find the first half in the dictionary, we're sure to
    find the second half

    :param part:
    :param key_length:
    :param index:
    :return:
    """

    # split reading into 2 (key_length)-char parts
    read1 = part[0:key_length]
    read2 = part[key_length:key_length + key_length]

    found = index.get(read1, 0)

    if found:
        return found
    else:
        return index.get(read2, 0) - key_length


def main():
    """
    Find a few short reads in a longer sequence,
    expecting at most one error per read
    :return:
    """

    key_length = 7
    reference = get_sequence()
    index = build_sequence_index(reference, key_length)

    assert get_position('ccggcctcgggaag', key_length, index) == 36
    assert get_position('ttgcggacgctagc', key_length, index) == 162
    assert get_position('tcgggctccccccg', key_length, index) == 87
    assert get_position('ggggggaaggcgga', key_length, index) == 114
    assert get_position('tctgtccccccccg', key_length, index) == 19

    print("All tests passed")


if __name__ == "__main__":
    main()
