from more_itertools import bucket, chunked, powerset, pairwise, filterfalse


def read_words(filename):
    with open(filename) as f:
        return {word.strip().lower() for word in f}


def create_graph(words):
    graph = {word: set() for word in words}
    for word in words:
        for i, char in product(range(len(word) + 1), string.ascii_lowercase):
            new_word = word[:i] + char + word[i:]
            if new_word in graph and new_word != word:
                graph[word].add(new_word)
                graph[new_word].add(word)
        for i, char in product(range(len(word)), string.ascii_lowercase):
            new_word = word[:i] + char + word[i+1:]
            if new_word in graph and new_word != word:
                graph[word].add(new_word)
                graph[new_word].add(word)
    return graph


def find_groups(graph):
    groups = []
    seen = set()
    for word in graph:
        if word not in seen:
            group = set(bucket(powerset(graph[word]), 2))
            seen.update(group)
            groups.append(group)
    return groups


def find_largest_group(groups):
    return max(groups, key=len)


def main():
    words = read_words('/usr/share/dict/words')
    graph = create_graph(words)
    groups = find_groups(graph)
    largest_group = find_largest_group(groups)
    print('Largest group:', largest_group)


if __name__ == '__main__':
    main()
