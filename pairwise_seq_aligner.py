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
        values = []

        # We create +1 column and +1 row because the 1x1 cell of the table
        # must contain 0 according to the algorithm
        for i in range(len(self.seq2)+1):
            # Create a row
            values.append([])
            for j in range(len(self.seq1)+1):
                # Create a column
                # If it is the first column, we assign gap penalties
                if i == 0:
                    values[i].append(j * self.gap)
                else:
                    # If it is the first row, we assign gap penalties
                    if j == 0:
                        values[i].append(i * self.gap)
                    else:
                        values[i].append(0)

        return values


    def _table_filling(self, values):
        """
        Fills the table using the Needleman-Wunsch recurrence relation.
        """
        rows, cols = len(values), len(values[0])

        # Create traceback direction matrix
        directions = [[""] * cols for _ in range(rows)]

        for i in range(1, len(self.seq2)+1):
            for j in range(1, len(self.seq1)+1):
                # Select max score and assign direction
                values[i][j], directions[i][j] = max(
                    (values[i - 1][j - 1] + self.matrix[self.seq2[i - 1], self.seq1[j - 1]], 'd'),  # Diagonal (match/mismatch)
                    (values[i - 1][j] + self.gap, 't'),     # Top (gap in seq1)
                    (values[i][j - 1] + self.gap, 'l')      # Left (gap in seq2)
                )

        return values, directions


    def _traceback(self, values, directions):
        """
        Traces back from the bottom-right of the DP table to reconstruct alignment.
        """
        aln1 = ""
        aln2 = ""

        i, j = len(self.seq2), len(self.seq1)  # Start from bottom-right

        while i > 0 or j > 0:
            if directions[i][j] == 'd':  # Diagonal (match/mismatch)
                aln1 = self.seq1[j - 1] + aln1
                aln2 = self.seq2[i - 1] + aln2
                i -= 1
                j -= 1
            elif directions[i][j] == 't':  # Up (gap in seq1)
                aln1 = "-" + aln1
                aln2 = self.seq2[i - 1] + aln2
                i -= 1
            else:  # Left (gap in seq2)
                aln1 = self.seq1[j - 1] + aln1
                aln2 = "-" + aln2
                j -= 1

        return aln1, aln2, values[-1][-1]


    def needleman_wunsch(self):
        """
        Traces back from the bottom-right of the DP table to reconstruct alignment.
        """
        values = self._initialization()
        values, directions = self._table_filling(values)
        return self._traceback(values, directions)

if __name__ == "__main__":
    aligner = PairwiseSeqAligner(
        "TGGGGGG",
        "TAGGTTAGG",
        substitution_matrix.SubstitutionMatrix("matrices/TTM.txt"),
        -2,
    )
    aln1, aln2, score = aligner.needleman_wunsch()
    print("Alignment 1: {0}, \nAlignment 2: {1}, \nScore: {2}".format(aln1, aln2, score))
