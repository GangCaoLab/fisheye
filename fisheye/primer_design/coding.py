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


class LLHC(HuffmanCoding):
    def __init__(self, *args, desire_depth=3, **kwargs):
        super().__init__(*args, **kwargs)
        self.desire_depth = desire_depth

    def reset_stack(self):
        stack = []
        def walk(node):
            if node.is_balance():
                stack.append(node)
            else:
                for c in node.childs:
                    walk(c)
        walk(self.tree.root)
        stack.sort(key=lambda n: (n.weight, n.value))
        self.tree.stack = stack

    def reduce_depth(self, diff):
        stack = self.tree.stack
        while diff > 0:
            n_childs = min(len(stack), len(self.code_book))
            childs = stack[:n_childs]
            stack = stack[n_childs:]
            new = Node.new_with_childs(childs)
            idx = bisect.bisect([(n.weight, n.value) for n in stack], (new.weight, new.value))
            stack.insert(idx, new)
            diff -= 1
        self.tree.stack = stack

    def coding(self, freq):
        hc = HuffmanCoding(self.code_book)
        self.tree = hc.coding(freq)
        if len(freq) > len(self.code_book) ** self.desire_depth:
            raise ValueError(f"Input can't coding with a {self.desire_depth} depth tree.")
        depth_diff = self.tree.root.depth - self.desire_depth
        if depth_diff > 0:
            self.reset_stack()
            self.reduce_depth(depth_diff)
            hc.tree = self.tree
            return hc.coding(freq, add_leaf=False)
        else:
            return self.tree


if __name__ == '__main__':
    import random
    freq = {}
    for i in range(1,1000) :
        freq[f'Gene{i}'] = random.randint(1,10)
    for i in range(1001,2100) :
        freq[f'Gene{i}'] = random.randint(100,1000)

    code_book = ['AA','AT','AG','AC','TT','TG','TC','GG','CG','CC']
    hc = HuffmanCoding(code_book)
    tree = hc.coding(freq)
    codes = hc.get_codes()
    c_codes = code_completion(codes)
    ##print(tree.root)
    ##print(tree.root.childs)

    ##llhc = LLHC(code_book, desire_depth=4)
    ##tree_ll = llhc.coding(freq)
    ##codes_ll = llhc.get_codes()
    ##c_codes_ll = code_completion(codes)
    ##print(tree_ll.root)
    ##print(tree_ll.root.childs)

    a = [1,1,1,1,1,1,1,1,1,1,5,5]
    freq = {f"g{i}": f for i, f in enumerate(a)}
    llhc = LLHC(['0', '1'], desire_depth=4)
    t = llhc.coding(freq)
