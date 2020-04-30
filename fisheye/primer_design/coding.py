import uuid
import bisect
from collections import namedtuple
import typing as t
from functools import reduce
import heapq


class Node(object):
    def __init__(self, id, childs=None, weight=0, depth=0, height=0):
        self.id = id
        self.weight = weight
        if childs is None:
            self.childs = []
        else:
            self.childs = childs
        self.depth = depth
        self.height = height

    def __repr__(self):
        return f"<Node id={self.id} val={self.weight} " +\
               f"depth={self.depth}, height={self.height}>"

    @staticmethod
    def new_with_childs(childs):
        id_ = str(uuid.uuid1())
        weight = 0; depths = []; heights = []
        for c in childs:
            weight += c.weight
            depths.append(c.depth)
            heights.append(c.height)
        depth = min(depths) - 1
        height = max(heights) + 1
        new = Node(id_, weight=weight, childs=childs,
                   depth=depth, height=height)
        return new


class Tree(object):
    def __init__(self):
        self.stack = []
        self.leafs = set()

    def add_leaf(self, node):
        self.stack.append(node)
        self.leafs.add(node)

    def pop_min(self, key=lambda n: n.weight):
        min_val = min(self.stack, key=key)
        min_idx = self.stack.index(min_val)
        return self.stack.pop(min_idx)

    def push(self, node):
        self.stack.append(node)

    @property
    def root(self):
        return self.stack[-1]


def code_completion(codes, code_length=6, code_unit_length=2):
    res = {}
    for gene, code in codes.items():
        c_len = len(code)
        if c_len < code_length:
            len_diff = code_length - c_len
            suffix = code[-code_unit_length:] * (len_diff//code_unit_length)
            res[gene] = code + suffix
        else:
            res[gene] = code
    return res


Item = namedtuple('Item', ['name', 'weight', 'd'])


class Package(object):

    def __init__(self, set_:t.Set[Item]):
        self.set = set_

    def merge(self, other: 'Package'):
        return Package(self.set | other.set)

    @property
    def weight(self):
        return sum([i.weight for i in self.set])


Number = t.Union[int, float]


class LLHC(object):
    """Length-limited huffman coding by Package-Merge algorithm.

    Reference:
    [1] Larmore, L. L. & Hirschberg, D. S. A fast algorithm for optimal
        length-limited Huffman codes. J. ACM 37, 464–473 (1990).
    """
    def __init__(self, code_book: t.List[str], L: int):
        self.code_book = code_book
        self.L = L
        self._tree = None

    def coding(self, freq: t.Mapping[str, Number]) -> t.Mapping[str, str]:
        if len(freq) > (len(self.code_book) ** self.L):
            raise ValueError(f"Input(size={len(freq)}) can't encode with a {self.L} depth " + \
                             f"{len(self.code_book)}-branch coding tree.")
        nodeset = self.get_nodeset(freq)
        self._nodeset = nodeset
        tree = self.nodeset2tree(nodeset)
        self._tree = tree
        codes = self.decode_tree(tree)
        return codes

    def get_nodeset(self, freq):
        L = self.L
        I = [(k, v) for k, v in freq.items()]
        P = [set()]
        for d in range(1, L+1):
            P.append({Package({Item(i[0], i[1], d)}) for i in I})

        for d in range(L, 0, -1):
            while len(P[d]) >= len(self.code_book):
                pkgs = []
                for i in range(len(self.code_book)):
                    p = min(P[d], key=lambda p: p.weight)
                    P[d].remove(p)
                    pkgs.append(p)
                p_new = reduce(lambda a,b: a.merge(b), pkgs)
                P[d-1].add(p_new)

        S = reduce(lambda a, b: a.merge(b), P[0]).set
        nodeset = sorted(S, key=lambda t: (t.name, t.d))
        return nodeset

    def nodeset2tree(self, nodeset: t.List[Item]):
        # only keep item with max d
        maxd_items = []
        old = None
        for i in nodeset:
            if old and (i.name != old.name):
                maxd_items.append(old)
            old = i
        maxd_items.append(i)
        maxd_items.sort(key=lambda i: i.d)

        # create tree, init leafs
        tree = Tree()
        for i in maxd_items:
            tree.add_leaf(Node(i.name, weight=i.weight, depth=i.d))
        
        # build tree
        while len(tree.stack) > 1:
            n_childs = min(len(tree.stack), len(self.code_book))
            childs = []
            for _ in range(n_childs):
                node = tree.pop_min(key=lambda n: (-n.depth, n.weight))
                childs.append(node)
            new = Node.new_with_childs(childs)
            tree.push(new)

        return tree

    def decode_tree(self, tree):
        codes = {}
        def walk(root, c=""):
            if not root.childs:  # leaf node
                codes[root.id] = c
            else:
                for i, child in enumerate(root.childs):
                    walk(child, c+self.code_book[i])
        walk(tree.root)
        return codes




if __name__ == '__main__':
    #freq = [100,1,1,1,1,1,1,1,1,1,1,1,5,5]
    #freq = [1, 1, 3, 4, 8, 100]
    #freq = {f"g{i}":f for i, f in enumerate(freq)}
    #llhc = LLHC(['0', '1'], 10)
    #codes = llhc.coding(freq)
    #print(codes)

    import random
    freq = {}
    for i in range(1,100) :
        freq[f'Gene{i}'] = random.randint(1,10)
    for i in range(101,210) :
        freq[f'Gene{i}'] = random.randint(100,1000)
    cb = ['AA','AT','AG','AC','TT','TG','TC','GG','CG','CC']
    llhc = LLHC(cb, 3)
    codes = llhc.coding(freq)
    print(codes)

