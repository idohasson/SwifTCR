import itertools
import collections

def generate_skipgrams(word):
    skipgrams = set()
    for i in range(len(word)):
        skipgrams.update(''.join(combo) for combo in itertools.combinations(word, i+1))
    return skipgrams

def generate_word_groups(words):
    word_groups = collections.defaultdict(list)
    for word in words:
        skipgrams = generate_skipgrams(word)
        for skipgram in skipgrams:
            word_groups[skipgram].append(word)
    return word_groups

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if not s2:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, ch1 in enumerate(s1):
        current_row = [i + 1]
        for j, ch2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (ch1 != ch2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def generate_distance_groups(words, distance_func):
    word_groups = generate_word_groups(words)
    distance_groups = set()
    for group in word_groups.values():
        if len(group) > 1:
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    if distance_func(group[i], group[j]) == 1:
                        distance_groups.add((group[i], group[j]))
    return distance_groups

words = ["cat", "bat", "rat", "car", "bar", "far", "hello", "world", "hi", "bye"]

# Compute the Hamming distance between all pairs of words
hamming_groups = generate_distance_groups(words, hamming_distance)

# Compute the Levenshtein distance between all pairs of words
levenshtein_groups = generate_distance_groups(words, levenshtein_distance)

print("Hamming groups:", hamming_groups)
print("Levenshtein groups:", levenshtein_groups)
