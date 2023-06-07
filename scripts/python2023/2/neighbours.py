import string

def neighbours(word):
    for j in range(len(word)):
        for d in string.ascii_lowercase:
            word1 = ''.join(d if i == j else c for i, c in enumerate(word))
            if word1 != word and word1 in words:
                yield word1

def graph2(words):
    return {(w1, w2) for w1 in words for w2 in neighbours(w1) if w1 < w2}

words = ["cat", "bat", "rat", "car", "bar", "far", "hello", "world", "hi", "bye"]

# Compute the Hamming distance between all pairs of words
hamming_pairs = graph2(words)

print("Hamming pairs:", hamming_pairs)
