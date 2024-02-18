# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing

        # initialize matrix private attributes

        gap_open_penalty = self.gap_open
        gap_extend_penalty = self.gap_extend
        n = len(self._seqA)
        m = len(self._seqB)

        # initialize matrices for alighment score, gap, and backtrace all with same dimension
    
        self._align_matrix = np.zeros((m + 1, n + 1))
        self._score_matrix = np.zeros((m + 1, n + 1))
        self._gapA_matrix = np.zeros((m + 1, n + 1))
        self._gapB_matrix = np.zeros((m + 1, n + 1))

        # initialize starting values including gap penalties

        for i in range(1, 1 + m):
            self._align_matrix[i][0] = -np.inf
            self._gapA_matrix[i][0] = gap_open_penalty + (i * gap_extend_penalty)
            self._gapB_matrix[i][0] = -np.inf

        for j in range(1, 1 + n):
            self._align_matrix[0][j] = -np.inf
            self._gapA_matrix[0][j] = -np.inf
            self._gapB_matrix[0][j] = gap_open_penalty + (j * gap_extend_penalty)

        # TODO: Implement global alignment here

        # construct alighment matrix
            
        for i in range(1, 1 + m):
            for j in range(1, 1 + n):
                self._align_matrix[i][j] = max(
                    self._align_matrix[i - 1][j - 1] + self.sub_dict[(seqA[j - 1], seqB[i - 1])],
                    self._gapA_matrix[i - 1][j - 1] + self.sub_dict[(seqA[j - 1], seqB[i - 1])],
                    self._gapB_matrix[i - 1][j - 1] + self.sub_dict[(seqA[j - 1], seqB[i - 1])]
                )

        # construct gapA matrix
                self._gapA_matrix[i][j] = max(
                    (self._align_matrix[i - 1][j] + gap_open_penalty + gap_extend_penalty),
                    (self._gapA_matrix[i - 1][j] + gap_extend_penalty)
                )
                
        # construct gapB matrix
                self._gapB_matrix[i][j] = max(
                    (self._align_matrix[i - 1][j] + gap_open_penalty + gap_extend_penalty),
                    (self._gapB_matrix[i - 1][j] + gap_extend_penalty)
                )
                
        # construct score matrix from max of alignment, gapa, amd gapb matrices
                self._score_matrix[i][j] = max(self._align_matrix[i][j], self._gapA_matrix[i][j], self._gapB_matrix[i][j])


        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # initialize a traceback matrix

        n = len(self._seqA)
        m = len(self._seqB)

        traceback_matrix = np.ones((m + 1, n + 1)) * -np.inf

        # iterate through traceback matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
        
        # 0 is a match
                if self._score_matrix[i][j] == self._align_matrix[i][j]:
                    traceback_matrix[i][j] = 0

        # -1 is a gap in A
                elif self._score_matrix[i][j] == self._gapA_matrix[i][j]:
                    traceback_matrix[i][j] = -1

        # 1 is a gap in B
                elif self._score_matrix[i][j] == self._gapB_matrix[i][j]:
                    traceback_matrix[i][j] = 1

        # if function breaks then raise error
                else:
                    raise ValueError("Iteration through traceback matrix has failed.")

        # trace the sequence through the matrix and produce aligned seqA and seqB
        while i > 0 and j > 0:
            trace = traceback_matrix[i][j]
        
        # add match to both sequences
            if trace == 0:
                self.seqA_align += self._seqA[j - 1]
                self.seqB_align += self._seqB[i - 1]
                i -= 1
                j -= 1

        # add gap in A
            elif trace == -1:
                self.seqA_align += '-'
                self.seqB_align += self._seqB[i - 1]
                i -= 1

        # add gap in B
            elif trace == 1:
                self.seqA_align += self._seqA[j - 1]
                self.seqB_align += '-'
                j -= 1
        
        # terminate traceback when sequence A or B ends
        if i > 0 and j == 0:
            self.seqA_align += '-'
            self.seqB_align += self._seqB[i - 1]
            i -= 1

        if i == 0 and j > 0:
            self.seqA_align += self._seqA[j - 1]
            self.seqB_align += '-'
            j -= 1

        # invert the traceback to original sequence orientation
        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1]

        # calculate the alignment score
        num_gaps = 0
        for i in range(min(len(self.seqA_align), len(self.seqB_align))):
        
        # if A and B match, add subdict alighment score and reset num_gaps
            if self.seqA_align[i] == self.seqB_align[i]:
                self.alignment_score += self.sub_dict[(self.seqA_align[i], self.seqB_align[i])]
                num_gaps = 0
        
        # if there is a gap, add open penalty AND extend penalty if the gap is starting
        # otherwise only add extend penalty
            else:
                if self.seqA_align[i] == '-' or self.seqB_align[i] == '-':
                    if num_gaps == 0:
                        self.alignment_score += self.gap_open + self.gap_extend
                        num_gaps += 1
                    else:
                        self.alignment_score += self.gap_extend

        # if there is a mismatch, add mismatch penalty and reset num_gaps
                else:
                    self.alignment_score += self.sub_dict[(self.seqA_align[i], self.seqB_align[i])]
                    num_gaps = 0

        # return self alignment score, seqA and seqB
        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
