class Node:
    def __init__(self, val):
        self.val = val
        self.edges = {}

    def __lshift__(self, node):
        self.edges[node.val] = True

    def connected(self, node):
        return self.edges.get(node.val, False)

    def __repr__(self):
        return f"Val: {self.val}, edges: {', '.join(self.edges.keys())}"

class Graph:
    def __init__(self):
        self.vertices = {}

    def __lshift__(self, val):
        self.vertices[val] = Node(val)

    def connect(self, node1, node2):
        node1 << node2
        node2 << node1

    def __iter__(self):
        for val, node in self.vertices.items():
            yield val, node

    def get(self, val):
        return self.vertices.get(val, None)

CHARACTERS = list('abcdefghijklmnopqrstuvwxyz')
graph = Graph()

# ~ 240 000 words
with open('/usr/share/dict/words', 'r') as f:
    for word in f:
        word = word.strip().lower()
        graph << word

for val, node in graph:
    for char in CHARACTERS:
        for i in range(len(val)+1):
            node2 = graph.get(val[:i] + char + val[i:])
            if node2:
                graph.connect(node, node2)
            if i < len(val):
                node2 = graph.get(val[:i] + char + val[i+1:])
                if node2:
                    graph.connect(node, node2)
