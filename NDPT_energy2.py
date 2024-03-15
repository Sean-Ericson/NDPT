from decimal import InvalidOperation
import numpy as np
import math
import pickle

# rotate list ls to the right by n (negative to go left)
def rotate_list(ls, n):
    return ls[-n:] + ls[:-n]

# Determine if two lists are rotated versions of each other
def list_cyclic_equal(l1, l2):
    # Normal list equal
    def list_equal(l1, l2):
        for i in range(len(l1)):
            if l1[i] != l2[i]:
                return False
        return True

    if len(l1) != len(l2):
        return False

    for i in range(len(l1)):
        if list_equal(l1, rotate_list(l2, i)):
            return True
    return False

# Determine if two lists are equal as sets. (can't use set() if the elements aren't hashable)
def list_set_equal(l1, l2):
    return np.all(np.array([x in l2 for x in l1])) and np.all(np.array(x in l1 for x in l2))

# Return all first/last-term-canceling partitions for correction n
def smart_partition(n):
    # returns all *lists* n non-negative numbers that sum to m
    def partition(n,m):
        res = []
        if m == 0:
            return [[0 for _ in range(n)]]
        if n == 1:
            return [[m]]
        for i in range(m+1):
            for ls in partition(n-1, m-i):
                res += [[i] + ls]
        return res

    res = []
    # zero on each side
    for ls in partition(n-1, n-1):
        res.append([0] + ls + [0])
    # nonzero on each side
    for i in range(1, n):
        for j in range(1, n-i):
            for ls in partition(n-1, n-1-i-j):
                res.append([i] + ls + [j])
    return res

class MultiSet:
    def __init__(self) -> None:
        self.items = []

    def add(self, new_item, count, auto_clear=True):
        need_to_append = True
        for i in range(len(self.items)):
            current_item, current_count = self.items[i]
            if current_item == new_item:
                self.items[i] = (current_item, current_count + count)
                need_to_append = False
                break
        if need_to_append:
            self.items.append(tuple((new_item, count)))
        if auto_clear:
            self.clear_zero_count_items()

    def add_many(self, items, counts):
        for item, count in zip(items, counts):
            self.add(item, count, auto_clear=False)
        self.clear_zero_count_items()

    def clear_zero_count_items(self):
        self.items = [i for i in self.items if i[1] != 0]

    def __eq__(self, other):
        return list_set_equal(self.items, other.items)

    def __str__(self):
        return " ".join(["{}^{:d}".format(i[0], i[1]) if i[1] != 1 else "{}".format(i[0]) for i in self.items])

    def __repr__(self):
        return str(self)

class SigmaTerm:
    def __init__(self, indices):
        self.indices = indices

    def __eq__(self, other: object) -> bool:
        return len(self.indices) == len(other.indices) and np.all(np.array(self.indices) == np.array(other.indices))
    
    def __str__(self):
        return "Sigma({})".format(self.indices)

    def __repr__(self) -> str:
        return str(self)

    def ToLatex(self):
        return "\\Sigma_{{}}".format(self.indices)

class PerturbativeTerm:
    def __init__(self, v_exp, sigmas) -> None:
        self.v_exp = v_exp
        self.sigmas = sigmas

    @staticmethod
    def FromPartition(partition):
        # Combine first/last
        partition = [partition[0] + partition[-1]] + partition[1:-1]

        #rotate to start at a boundary 0
        for i in range(len(partition)):
            if partition[i] != 0 and partition[(i+1)%len(partition)] == 0:
                shift = i+1    
        partition = rotate_list(partition, -shift)
        
        pterm = PerturbativeTerm(0, MultiSet())
        
        # get info for v_exp
        zero_run_lens = [1]
        on_run = True
        for i in range(1, len(partition)):
            if partition[i] == 0:
                if on_run:
                    zero_run_lens[-1] += 1
                else:
                    zero_run_lens.append(1)
                    on_run = True
            else:
                on_run = False
        pterm.v_exp = sum(zero_run_lens) - len(zero_run_lens)

        # get sigmas
        number_groups = []
        in_group = False
        for i in range(1, len(partition)):
            if partition[i] != 0:
                if in_group:
                    number_groups[-1].append(partition[i])
                else:
                    number_groups.append([partition[i]])
                    in_group = True
            else:
                in_group = False
        pterm.sigmas.add_many(*zip(*[(SigmaTerm(group), 1) for group in number_groups]))
        
        return pterm

    def __eq__(self, other):
        return self.v_exp == other.v_exp and self.sigmas == other.sigmas

    def __str__(self) -> str:
        return "(V_00 ^ {:d}) ".format(self.v_exp) + str(self.sigmas) if self.v_exp > 0 else ""
    
    def string_with_coeff(self, coeff):
        return "{:d} * (V_00 ^ {:d}) ".format(coeff, self.v_exp) + str(self.sigmas) if self.v_exp > 0 else "{:d} * ".format(coeff) + str(self.sigmas)

    def __repr__(self) -> str:
        return str(self)

class EnergyCorrection:
    def __init__(self, n):
        self.n = n
        self.SigmasByV00 = dict()

    def calc(self):
        # Generate perms and combine first/last 
        p_terms = MultiSet()
        foo = smart_partition(self.n)
        bar = "eyyyylmao"
        print(len(foo))
        for part in foo:
            p_terms.add(PerturbativeTerm.FromPartition(part), (-1)**(part.count(0)), auto_clear=False)
        p_terms.clear_zero_count_items()

        #Sort by V00
        self.v_max = max(list(zip(*p_terms.items))[0], key=lambda t: t.v_exp).v_exp
        for i in range(self.v_max + 1):
            self.SigmasByV00[i] = [t for t in p_terms.items if t[0].v_exp == i]

        self.p_terms = p_terms

    def print_pterms_by_v00(self):
        for i in range(self.v_max + 1):
            for t in [t for t in self.p_terms.items if t[0].v_exp == i]:
                print(t[0].string_with_coeff(t[1]))
        
    def to_list(self):
        return [(v_exp, [([(sigma.indices, exp) for sigma,exp in sigmas.sigmas.items], coeff) for sigmas,coeff in pterm]) for v_exp, pterm in self.SigmasByV00.items()]
    
for i in range(2,13):
    correction = EnergyCorrection(i)
    correction.calc()
    data = correction.to_list()
    print("Terms in correction {}: {}".format(i, len(correction.p_terms.items)))
    #with open("Coded_Corrections.pkl", 'wb') as file:
    #    pickle.dump((i, data), file)