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
        Fills the table using the recurrence relation.
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
        Needleman-Wunsch algorithm for global alignment:

        1. Initialization
            In this step a table is produced.
            The rows of the table represent characters of the second sequence.
            The columns of the table represent characters of the first sequence.
            The 0x0 element is always a zero because according to the algorithm
            the first character of the second sequence is represented by the second row
            and first character of the first sequence is represented by the second column.
            Also, the gap penalties are assigned for all elements of the first row and all elements of the first column.
        2. Table filling
            In this step the table is filled with values and directions for tracing back.
            We iterate over elements and assign a value by choosing a maximum of following options:
                - The value of previous diagonal element + match/mismatch score (taken from a matrix)
                - The value of the element from the left + gap penalty
                - The value of the element from the top + gap penalty
            The last element is the maximum value
        3. Traceback
            In this step we collect the alignments by tracing back the table.
            We start from the last element and use directions to go over elements to the first elements:
                - If a value has the diagonal direction, we know that we should go to the previous diagonal element and add to our alignments
                    the characters corresponding to the previous diagonal element.
                - If a value has the left direction, we go to the left, add to the first alignment the character corresponding
                    to the previous left element and add to the second alignment the gap character.
                - If a value has the top direction, we go to the top, add to the second alignment the character corresponding
                    to the previous left element and add to the first alignment the gap character.
        """
        modifier = lambda a: a
        table = self._initialization(modifier)
        table = self._table_filling(table, modifier, False)
        condition = lambda t,a,b: True if a > 0 or b > 0 else False
        return self._traceback(table, len(self.seq2), len(self.seq1), table[-1][-1][0], condition)

    def smith_waterman(self):
        """
        Smith-Waterman algorithm for local alignment:

        1. Initialization
            In this step a table is produced.
            The rows of the table represent characters of the second sequence.
            The columns of the table represent characters of the first sequence.
            The 0x0 element is always a zero because according to the algorithm
            the first character of the second sequence is represented by the second row
            and first character of the first sequence is represented by the second column.
            The gap penalties are not assigned because according to the algorithm the minimum score is zero,
            so, we just assign zero to all elements of the first row and all elements of the first column.
        2. Table filling
            In this step the table is filled with values and directions for tracing back.
            We iterate over elements and assign a value by choosing a maximum of following options:
                - The value of previous diagonal element + match/mismatch score (taken from a matrix)
                - The value of the element from the left + gap penalty
                - The value of the element from the top + gap penalty
            In case of parity in we should override direction to '' that means that local alignment is over here.
            The last element is the maximum value
        3. Traceback
            In this step we collect the alignments by tracing back the table.
            First of all we iterate over the table to find the element with the maximum score.
            We start from this and use directions to go over elements until we encounter zero direction (''):
                - If a value has the diagonal direction, we know that we should go to the previous diagonal element and add to our alignments
                    the characters corresponding to the previous diagonal element.
                - If a value has the left direction, we go to the left, add to the first alignment the character corresponding
                    to the previous left element and add to the second alignment the gap character.
                - If a value has the top direction, we go to the top, add to the second alignment the character corresponding
                    to the previous left element and add to the first alignment the gap character.

        Differences with NW algorithm:
            - The positive constraint: all values must be >= 0
            - A direction might be '' meaning that local alignment is over here
            - We start collect alignments from the element with the maximum score and finish when we encounter zero direction ('')
            - Output is not required to contain whole sequences, it contains the best local alignments
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
