# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    score, s1_s2_align, s2_s1_align = nw.align(seq1, seq2)

    assert s1_s2_align == 'MYQR'

    assert s2_s1_align == 'M-QR'
    
    print(nw._align_matrix)
    print(nw._gapA_matrix)
    print(nw._gapB_matrix)
    print(type(nw._align_matrix[0][0]))
    assert np.all(nw._align_matrix[0] == [0, -np.inf, -np.inf, -np.inf, -np.inf])

    assert np.all(nw._gapA_matrix == [
            [-10, -np.inf, -np.inf, -np.inf],
            [-11, -12, -6, -7],
            [-12, -13, -14, -7],
            [-13, -14, -15, -12],
            [-14, -15, -16, -17],
        ]
    )

    assert np.all(nw._gapB_matrix == [
            [-10, -11, -12, -13],
            [-np.inf, -12, -13, -14],
            [-np.inf, -6, -14, -15],
            [-np.inf, -7, -7, -16],
            [-np.inf, -8, -8, -6],
        ]
    )
def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    score, s3_s4_align, s4_s3_align = nw.align(seq3, seq4)
    
    assert score == 17
    assert s3_s4_align == 'MAVHQLIRRP'
    assert s4_s3_align == 'M---QLIRHP'



