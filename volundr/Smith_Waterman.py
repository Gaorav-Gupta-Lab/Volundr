# -*- coding: utf-8 -*-
"""
Smith_Waterman.py v0.1.0
    Dec. 1, 2016
    Dennis A. Simpson
    Changed versioning to conform to semantic versioning (http://semver.org/).
Smith_Waterman.py v1.0
    Dennis A. Simpson
    August 19, 2016
    A Python implementation of the Smith-Waterman algorithm for local alignment of nucleotide sequences.Based on (c) 2013
    Ryan Boehning https://gist.github.com/radaniba/11019717.  Modified to be a class for use in Volundr.  The code still
    needs refactoring.  A decision needs to be made as to keeping the string alignment outputs.
"""

__author__ = 'Dennis A. Simpson'
__version__ = '0.1.0'
__package__ = 'Volundr'


class SmithWaterman:
    def __init__(self):
        """
        The values for match, mismatch, and gap are from Wikipedia.
        en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        """
        self.match = 2
        self.mismatch = -1
        self.gap = -1
        self.seq1 = None
        self.seq2 = None
        self.max_pos = None
        self.score_matrix = None
        self.aligned_seq1 = []
        self.aligned_seq2 = []

    def get_score(self, seq1, seq2):

        self.seq1 = seq1
        self.seq2 = seq2
        SmithWaterman.create_score_matrix(self)
        SmithWaterman.traceback(self)

        # Build the string as a list of characters to avoid costly string
        # concatenation.
        idents, gaps, mismatches = 0, 0, 0
        algnment_string = []

        for base1, base2 in zip(self.aligned_seq1, self.aligned_seq2):
            if base1 == base2:
                algnment_string.append('|')
                idents += 1
            elif '-' in (base1, base2):
                algnment_string.append(' ')
                gaps += 1
            else:
                algnment_string.append(':')
                mismatches += 1

        # return ''.join(alignment_string), idents, gaps, mismatches
        return len(self.aligned_seq1), gaps+mismatches

    def create_score_matrix(self):
        """
        Create a matrix of scores representing trial alignments of the two sequences.
        Sequence alignment can be treated as a graph search problem. This function
        creates a graph (2D matrix) of scores, which are based on trial alignments
        of different base pairs. The path with the highest cumulative score is the
        best alignment.
        """

        # The scoring matrix contains an extra row and column for the gap (-), hence the +1 here.
        rows = len(self.seq1) + 1
        cols = len(self.seq2) + 1

        self.score_matrix = [[0 for col in range(cols)] for row in range(rows)]

        # Fill the scoring matrix.
        max_score = 0
        self.max_pos = None    # The row and column of the highest score in matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                score = SmithWaterman.__calc_score(self, i, j)
                if score > max_score:
                    max_score = score
                    self.max_pos = (i, j)

                self.score_matrix[i][j] = score

        assert self.max_pos is not None, 'the x, y position with the highest score was not found'

        return

    def __calc_score(self, x, y):
        """
        Calculate score for a given x, y position in the scoring matrix.
        The score is based on the up, left, and upper-left neighbors.
        """
        similarity = self.match if self.seq1[x - 1] == self.seq2[y - 1] else self.mismatch

        diag_score = self.score_matrix[x - 1][y - 1] + similarity
        up_score = self.score_matrix[x - 1][y] + self.gap
        left_score = self.score_matrix[x][y - 1] + self.gap

        return max(0, diag_score, up_score, left_score)

    def traceback(self):
        """
        Find the optimal path through the matrix.
        This function traces a path from the bottom-right to the top-left corner of
        the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
        or both of the sequences being aligned. Moves are determined by the score of
        three adjacent squares: the upper square, the left square, and the diagonal
        upper-left square.
        WHAT EACH MOVE REPRESENTS
            diagonal: match/mismatch
            up:       gap in sequence 1
            left:     gap in sequence 2
        """

        END, DIAG, UP, LEFT = range(4)

        x, y = self.max_pos  # start_pos
        move = next_move(self.score_matrix, x, y)
        while move != END:
            if move == DIAG:
                self.aligned_seq1.append(self.seq1[x - 1])
                self.aligned_seq2.append(self.seq2[y - 1])
                x -= 1
                y -= 1
            elif move == UP:
                self.aligned_seq1.append(self.seq1[x - 1])
                self.aligned_seq2.append('-')
                x -= 1
            else:
                self.aligned_seq1.append('-')
                self.aligned_seq2.append(self.seq2[y - 1])
                y -= 1

            move = next_move(self.score_matrix, x, y)

        self.aligned_seq1.append(self.seq1[x - 1])
        self.aligned_seq2.append(self.seq2[y - 1])
        ''.join(reversed(self.aligned_seq1))
        ''.join(reversed(self.aligned_seq2))

        assert len(self.aligned_seq1) == len(self.aligned_seq2), 'aligned strings are not the same size'

        # return ''.join(reversed(self.aligned_seq1)), ''.join(reversed(self.aligned_seq2))

# def main():
#
#     # Pretty print the results. The printing follows the format of BLAST results
#     # as closely as possible.
#     alignment_str, idents, gaps, mismatches = alignment_string(seq1_aligned, seq2_aligned)
#     alength = len(seq1_aligned)
#     print()
#     print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
#           alength, idents / alength, gaps, alength, gaps / alength))
#     print()
#     for i in range(0, alength, 60):
#         seq1_slice = seq1_aligned[i:i+60]
#         print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
#         print('             {0}'.format(alignment_str[i:i+60]))
#         seq2_slice = seq2_aligned[i:i+60]
#         print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
#         print()


def next_move(score_matrix, x, y):
    diag = score_matrix[x - 1][y - 1]
    up = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')


def alignment_string(aligned_seq1, aligned_seq2):
    """
    Construct a special string showing identities, gaps, and mismatches.
    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.
    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    """
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    algnment_string = []
    for base1, base2 in zip(aligned_seq1, aligned_seq2):
        if base1 == base2:
            algnment_string.append('|')
            idents += 1
        elif '-' in (base1, base2):
            algnment_string.append(' ')
            gaps += 1
        else:
            algnment_string.append(':')
            mismatches += 1

    return ''.join(algnment_string), idents, gaps, mismatches


def print_matrix(matrix):
    """
    Print the scoring matrix.
    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    """
    for row in matrix:
        for col in row:
            print('{0:>4}'.format(col))
        print()

"""
class ScoreMatrixTest(unittest.TestCase):
    '''Compare the matrix produced by create_score_matrix() with a known matrix.'''
    def test_matrix(self):
        # From Wikipedia (en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
        #                -   A   C   A   C   A   C   T   A
        known_matrix = [[0,  0,  0,  0,  0,  0,  0,  0,  0],  # -
                        [0,  2,  1,  2,  1,  2,  1,  0,  2],  # A
                        [0,  1,  1,  1,  1,  1,  1,  0,  1],  # G
                        [0,  0,  3,  2,  3,  2,  3,  2,  1],  # C
                        [0,  2,  2,  5,  4,  5,  4,  3,  4],  # A
                        [0,  1,  4,  4,  7,  6,  7,  6,  5],  # C
                        [0,  2,  3,  6,  6,  9,  8,  7,  8],  # A
                        [0,  1,  4,  5,  8,  8, 11, 10,  9],  # C
                        [0,  2,  3,  6,  7, 10, 10, 10, 12]]  # A

        global seq1, seq2
        seq1 = 'AGCACACA'
        seq2 = 'ACACACTA'
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        matrix_to_test, max_pos = create_score_matrix(rows, cols)
        self.assertEqual(known_matrix, matrix_to_test)
"""
