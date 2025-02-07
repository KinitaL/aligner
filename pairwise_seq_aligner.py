import substitution_matrix

class PairwiseSeqAligner:
    """
    PairSeqAligner is a pairwise sequence aligner that implement
        - Needleman-Wunsch algorithm for global alignment
        - Smith-Waterman algorithm for local alignment
    """
    def __init__(self, seq1, seq2, matrix, gap):
        self.seq1 = seq1
        self.seq2 = seq2
        self.matrix = matrix
        self.gap = gap


    def _initialization(self):
        """
        Initializes the matrix with gap penalties.
        """
        table = []

        # We create +1 column and +1 row because the 1x1 cell of the table
        # must contain 0 according to the algorithm
        for i in range(len(self.seq2)+1):
            # Create a row
            table.append([])
            for j in range(len(self.seq1)+1):
                # Create a column
                # If it is the first column, we assign gap penalties
                if i == 0:
                    table[i].append((j * self.gap, ''))
                else:
                    # If it is the first row, we assign gap penalties
                    if j == 0:
                        table[i].append((i * self.gap, ''))
                    else:
                        table[i].append((0, ''))

        return table


    def _table_filling(self, table):
        """
        Fills the table using the Needleman-Wunsch recurrence relation.
        """
        for i in range(1, len(self.seq2)+1):
            for j in range(1, len(self.seq1)+1):
                # Select max score and assign direction
                table[i][j] = max(
                    (table[i - 1][j - 1][0] + self.matrix[self.seq2[i - 1], self.seq1[j - 1]], 'd'),  # Diagonal (match/mismatch)
                    (table[i - 1][j][0] + self.gap, 't'),     # Top (gap in seq1)
                    (table[i][j - 1][0] + self.gap, 'l')      # Left (gap in seq2)
                )

        return table


    def _traceback(self, table):
        """
        Traces back from the bottom-right of the DP table to reconstruct alignment.
        """
        aln1 = ""
        aln2 = ""

        i, j = len(self.seq2), len(self.seq1)  # Start from bottom-right

        while i > 0 or j > 0:
            if table[i][j][-1] == 'd':  # Diagonal (match/mismatch)
                aln1 = self.seq1[j - 1] + aln1
                aln2 = self.seq2[i - 1] + aln2
                i -= 1
                j -= 1
            elif table[i][j][-1] == 't':  # Up (gap in seq1)
                aln1 = "-" + aln1
                aln2 = self.seq2[i - 1] + aln2
                i -= 1
            else:  # Left (gap in seq2)
                aln1 = self.seq1[j - 1] + aln1
                aln2 = "-" + aln2
                j -= 1

        return aln1, aln2, table[-1][-1][0]


    def needleman_wunsch(self):
        """
        Needleman-Wunsch algorithm for global alignment
        """
        table = self._initialization()
        table = self._table_filling(table)
        return self._traceback(table)

if __name__ == "__main__":
    aligner = PairwiseSeqAligner(
        "ATGTCCCGTAATA",
        "ATGTCGTAATATGC",
        substitution_matrix.SubstitutionMatrix("matrices/TTM.txt"),
        -2,
    )
    aln1, aln2, score = aligner.needleman_wunsch()
    print("Alignment 1: {0}, \nAlignment 2: {1}, \nScore: {2}".format(aln1, aln2, score))
