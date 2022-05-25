import random


char_to_index_map = {
    "A": 0, "a": 0,
    "C": 1, "c": 1,
    "G": 2, "g": 2,
    "T": 3, "t": 3
}

def modify_sequence(sequence, num_errors, substitution_percentage=0.95):
    """
    Randomly modify a sequence: add several mismatches or indels, according
    to the specified `num_errors`.

    Args:
        sequence: A DNA sequence.
        num_errors: The total number of errors to be introduced.
        substitution_percentage: the percentage of mismatches in all the errors
            (mismatch / (mismatch + indel)).
    """
    sequence_copy = sequence
    print("--", sequence_copy)
    for _ in range(num_errors):
        random_number = random.uniform(0, 1)
        if random_number < substitution_percentage:
            # add a mismatch in the sequence
            i = random.randint(0, len(sequence_copy) - 1)
            current_char = sequence_copy[i]
            candidate_char = ['A', 'C', 'G', 'T']
            candidate_char.remove(current_char)
            sequence_copy = sequence_copy[:i] + random.choice(candidate_char) + sequence_copy[i+1:]
        elif random_number < substitution_percentage + (1-substitution_percentage)/2:
            # add an insert
            i = random.randint(0, len(sequence_copy))
            sequence_copy = sequence_copy[:i] + random.choice(['A', 'C', 'G', 'T']) + sequence_copy[i:]
        else:
            # add a delete
            i = random.randint(0, len(sequence_copy) - 1)
            sequence_copy = sequence_copy[:i] + sequence_copy[i+1:]
        
        print("->", sequence_copy)
        
    return sequence_copy


if __name__ == "__main__":
    print(modify_sequence("ACGCT", 2))