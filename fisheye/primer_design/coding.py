class Node(object):
    def __init__(self, id, childs=None, value=0, depth=0):
        self.id = id
        self.value = value
        if childs is None:
            self.childs = []
        else:
            self.childs = childs
        self.depth = depth

    def __repr__(self):
        return f"""<Node id={self.id} val={self.value}>"""

    @staticmethod
    def new_with_childs(childs):
        id_ = "_".join([c.id for c in childs])
        value = sum([c.value for c in childs])
        depth = max([n.depth for n in childs]) + 1
        new = Node(id_, value=value, childs=childs, depth=depth)
        return new


class Tree(object):
    def __init__(self):
        self.stack = []
        self.leafs = set()
        self.nodes = set()

    def add_leaf(self, node):
        self.stack.append(node)
        self.leafs.add(node)
        self.nodes.add(node)

    def pop_min(self):
        min_val = min(self.stack, key=lambda n: n.value)
        min_idx = self.stack.index(min_val)
        return self.stack.pop(min_idx)

    def push(self, node):
        self.nodes.add(node)
        self.stack.append(node)

    @property
    def root(self):
        return self.stack[-1]

    @property
    def nonleafs(self):
        return {n for n in self.nodes if n.childs}




class HuffmanCoding(object):
    def __init__(self, code_book):
        self.code_book = code_book
        self.tree = Tree()

    def coding(self, freq):
        tree = self.tree
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


if __name__ == '__main__':
    import random
    freq = {}
    for i in range(1,100) :
        freq[f'Gene{i}'] = random.randint(1,10)
    for i in range(101,210) :
        freq[f'Gene{i}'] = random.randint(100,1000)

    code_book = ['AA','AT','AG','AC','TT','TG','TC','GG','CG','CC']
    hc = HuffmanCoding(code_book)
    tree = hc.coding(freq)
    codes = hc.get_codes()
    c_codes = code_completion(codes)
