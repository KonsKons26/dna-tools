
class DNAAssembler:
    def __init__(self, reads):
        self.reads = reads

    # given two sequences and an offset, count the number of matching bases
    def score(self, sequence1, sequence2, offset):
        return sum(
            1 for position in range(max(0 - offset, 0),
                                    min(len(sequence2) - offset, len(sequence2), len(sequence1) - offset))
            if sequence2[position] == sequence1[position + offset]
        )

    # given two sequences, find the offset which gives the best score
    def find_best_offset(self, sequence1, sequence2):
        return max(
            (self.score(sequence1, sequence2, offset), offset, sequence2, sequence1)
            for offset in range(1 - len(sequence2), len(sequence1))
        )

    # given a single sequence and a collection of others, find the other sequence with the best match score
    def find_best_match(self, sequence, others):
        return max(
            self.find_best_offset(sequence, sequence2)
            for sequence2 in others if sequence2 != sequence
        )

    # given two sequences and an offset, calculate the consensus
    def consensus(self, score, offset, sequence1, sequence2):
        return sequence2[0:max(0, offset)] + sequence1 + sequence2[len(sequence1) + offset:]

    # given a sequence and collection of others, return the complete consensus using recursion
    def assemble(self, sequence, others):
        if len(others) == 1:
            return self.consensus(*self.find_best_match(sequence, others))
        else:
            best_match = self.find_best_match(sequence, others)
            return self.assemble(
                self.consensus(*best_match),
                [y for y in others if y != best_match[2]]
            )

    # given a collection of sequences, call assemble() to start the recursion
    def assemble_helper(self):
        return self.assemble(self.reads[0], self.reads[1:])
