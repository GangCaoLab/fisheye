import uuid
import bisect


class Node(object):
    def __init__(self, id, childs=None, value=0, depth=0, weight=1):
        self.id = id
        self.value = value
        if childs is None:
            self.childs = []
        else:
            self.childs = childs
        self.depth = depth
        self.weight = weight

    def __repr__(self):
        return f"""<Node id={self.id} val={self.value} depth={self.depth} weight={self.weight}>"""

    @staticmethod
    def new_with_childs(childs):
        id_ = str(uuid.uuid1())
        value = 0; depths = []; weight = 1
        for c in childs:
            value += c.value
            weight += c.weight
            depths.append(c.depth)
        depth = max(depths) + 1
        new = Node(id_, value=value, childs=childs, depth=depth, weight=weight)
        return new

    def is_balance(self):
        b = True
        childs = self.childs
        for i in range(1, len(self.childs)):
            b &= (childs[i].weight == childs[i-1].weight)
            if not b:
                break
        return b


class Tree(object):
    def __init__(self):
        self.stack = []
        self.leafs = set()

    def add_leaf(self, node):
        self.stack.append(node)
        self.leafs.add(node)

    def pop_min(self):
        min_val = min(self.stack, key=lambda n: n.value)
        min_idx = self.stack.index(min_val)
        return self.stack.pop(min_idx)

    def push(self, node):
        self.stack.append(node)

    @property
    def root(self):
        return self.stack[-1]


class HuffmanCoding(object):
    def __init__(self, code_book):
        self.code_book = code_book
        self.tree = Tree()

    def coding(self, freq, add_leaf=True):
        tree = self.tree
        if add_leaf:
            for gene, fq in freq.items():
                tree.add_leaf(Node(gene, value=fq))

        while len(tree.stack) > 1:
            childs = []
            n_childs = min(len(tree.stack), len(self.code_book))
            for i in range(n_childs):
                node = tree.pop_min()
                childs.append(node)
            new = Node.new_with_childs(childs)
            tree.push(new)
        return tree

    def get_codes(self):
        codes = {}
        def walk(root, c=""):
            if not root.childs:  # leaf node
                codes[root.id] = c
            else:
                for i, child in enumerate(root.childs):
                    walk(child, c+self.code_book[i])
        walk(self.tree.root)
        return codes


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


def LLHC(freq, lmax):
    I = sorted([(k, v) for k, v in freq.items()],
               key=lambda t: t[1])

    def mincost(i, X, L):
        if (X >= 2**(-lmax)) and (i <= len(I)-1):
            vals = []; lens = []
            for j in range(1, lmax+1):
                c, l = mincost(i+1, X-2**(-j), L+(j,))
                vals.append(2**j * I[i][1] + c)
                lens.append(l)
            minval = min(vals)
            minidx = vals.index(minval)
            return minval, lens[minidx]
        elif (X == 0) and (i == len(I)):
            return 0, L
        else:
            return float('inf'), tuple()

    return I, mincost(0, 1, tuple())



if __name__ == '__main__':
    freq = [100,1,1,1,1,1,1,1,1,1,1,1,5,5]
    freq = {f"g{i}":f for i, f in enumerate(freq)}
    minc = LLHC(freq,4)
    print(minc)
    #import random
    #freq = {}
    #for i in range(1,1000) :
    #    freq[f'Gene{i}'] = random.randint(1,10)
    #for i in range(1001,2100) :
    #    freq[f'Gene{i}'] = random.randint(100,1000)

    #code_book = ['AA','AT','AG','AC','TT','TG','TC','GG','CG','CC']
    #hc = HuffmanCoding(code_book)
    #tree = hc.coding(freq)
    #codes = hc.get_codes()
    #c_codes = code_completion(codes)
    ##print(tree.root)
    ##print(tree.root.childs)

