import re
import parasail


op_consumes = {
        "M": (True, True),
        "=": (True, True),
        "I": (True, False),
        "D": (False, True),
        "N": (False, True),
        "S": (True, False),
        "H": (False, False),
        "P": (False, False),
        "X": (True, True)
    }

cigar_pattern = re.compile("({})".format("|".join(list(op_consumes.keys()))))


def get_and_prepare_cigar(seq1, seq2, open=11, extend=1, matrix=parasail.blosum62):

    """This function is a wrapper around the alignment functions from parasail.
    It further prepares the cigar, by turning it into a list of tuples (operation, length)

    """

    result = parasail.sg_trace_scan_sat(seq1, seq2, open, extend, matrix)
    values = re.split(cigar_pattern, result.cigar.decode.decode())
    values = [(int(values[_ * 2]), values[_ * 2 + 1]) for _ in range(int((len(values) - 1) / 2))]
    return result, values


def cigar_length_in_common(cigar):

    return sum(length for length, op in cigar if any(op_consumes[op]))


def create_translation_array(cigar):

    common_length = cigar_length_in_common(cigar)
    # print(common_length)
    query_array = [None] * common_length
    target_array = [None] * common_length

    # Now we have to translate the exons into this common space
    # We will do this by creating "exons" derived from the alignment

    common_pos = 0
    query_pos = 0
    target_pos = 0

    for length, op in cigar:
        consumes_query, consumes_target = op_consumes[op]
        if not any((consumes_query, consumes_target)):
            continue
        else:
            for pos in range(common_pos, common_pos + length):
                if consumes_query:
                    query_pos += 1
                    query_array[pos] = query_pos
                if consumes_target:
                    target_pos += 1
                    target_array[pos] = target_pos
            common_pos += length

    return query_array, target_array


def transfer_pos_and_exons(cigar, query_exons, target_exons):

    query_array, target_array = create_translation_array(cigar)

    for pos in range(1, query_exons[-1][1] + 1):
        assert pos in query_array, (
        pos, min([_ for _ in query_array if _ is not None]), max([_ for _ in query_array if _ is not None]))
    for pos in range(1, target_exons[-1][1] + 1):
        assert pos in target_array


def transfer_exons(exons, c_array):
    """This static method will translate one set of exons into a different coordinate system
    (given an array that translates each position in system 1 to a position in system 2)"""

    c_exons = []
    for exon in exons:
        start, end = exon
        c_start = c_array.index(start)
        c_end = c_array.index(end)
        c_exons.append((c_start, c_end))

    return c_exons