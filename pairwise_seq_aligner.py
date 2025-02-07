import substitution_matrix

class PairwiseSeqAligner:
    # TODO: add description for both algorithm and stress main differences between them (aim, workflow, output etc.)
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


    def _initialization(self, modifier):
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
                    table[i].append((modifier(j * self.gap), ''))
                else:
                    # If it is the first row, we assign gap penalties
                    if j == 0:
                        table[i].append((modifier(i * self.gap), ''))
                    else:
                        table[i].append((modifier(0), ''))

        return table


    def _table_filling(self, table, modifier, is_sw):
        """
        Fills the table using the Needleman-Wunsch recurrence relation.
        """
        for i in range(1, len(self.seq2)+1):
            for j in range(1, len(self.seq1)+1):
                # Select max score and assign direction
                value = max(
                    (modifier(table[i - 1][j - 1][0] + self.matrix[self.seq2[i - 1], self.seq1[j - 1]]), 'd'),  # Diagonal (match/mismatch)
                    (modifier(table[i - 1][j][0] + self.gap), 't'),     # Top (gap in seq1)
                    (modifier(table[i][j - 1][0] + self.gap), 'l')     # Left (gap in seq2)
                )
                # In case of parity in local alignment we should override direction to ''
                if is_sw and value[0] == 0:
                    value = (0, '')
                table[i][j] = value


        return table


    def _traceback(self, table, i, j, max_score, condition):
        """
        Traces back from the bottom-right of the DP table to reconstruct alignment.
        """
        aln1 = ""
        aln2 = ""

        while condition(table, i, j):
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

        return aln1, aln2, max_score


    def needleman_wunsch(self):
        """
        Needleman-Wunsch algorithm for global alignment
        """
        modifier = lambda a: a
        table = self._initialization(modifier)
        table = self._table_filling(table, modifier, False)
        condition = lambda t,a,b: True if a > 0 or b > 0 else False
        return self._traceback(table, len(self.seq2), len(self.seq1), table[-1][-1][0], condition)

    def smith_waterman(self):
        """
        Smith-Waterman algorithm for local alignment
        """
        modifier = lambda a: 0 if a<0 else a
        table = self._initialization(modifier)
        table = self._table_filling(table, modifier, True)
        max_score = -1
        for i in range(len(self.seq2) + 1):
            for j in range(len(self.seq1) + 1):
                if table[i][j][0] > max_score:
                    max_score = table[i][j][0]
                    start_i = i
                    start_j = j

        i = start_i
        j = start_j
        condition = lambda t,a,b: True if t[a][b] != '' and a > 0 and b > 0 else False
        return self._traceback(table, i, j, max_score, condition)

if __name__ == "__main__":
    aligner = PairwiseSeqAligner(
        "AGCTAAGCCTAGC",
        "AGCTGCGCAGCGAGCCTAGC",
        substitution_matrix.SubstitutionMatrix("matrices/TTM.txt"),
        -2,
    )
    aln1, aln2, score = aligner.needleman_wunsch()
    print("GLOBAL\nAlignment 1: {0}, \nAlignment 2: {1}, \nScore: {2}".format(aln1, aln2, score))
    aln1, aln2, score = aligner.smith_waterman()
    print("LOCAL\nAlignment 1: {0}, \nAlignment 2: {1}, \nScore: {2}".format(aln1, aln2, score))
