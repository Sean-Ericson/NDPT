from decimal import InvalidOperation
import numpy as np
import math

def cyclic_list_equal(l1, l2):
    if len(l1) != len(l2):
        return False
    n = len(l1)
    ks = [i for i in range(n) if l1[0] == l2[i]]
    if len(ks) == 0:
        return False
    for k in ks:
        equal = True
        for i in range(n):
            if l1[i] != l2[(k + i) % n]:
                equal = False
                break
        if equal:
            return True 
    return False

def set_lists_equal(l1, l2):
    def isSubset(l1, l2):
        for x in l1:
            if not (x in l2):
                return False
        return True
    return isSubset(l1, l2) and isSubset(l2, l1)

def insert(i, x, ls):
    l = list(ls)
    l.insert(i,x)
    return l

# return a list of lists, each of wich is ls w/ x inserted into a different place
def insert_into_every_pos(x, ls):
    return [insert(i,x, ls) for i in range(len(ls)+1)]

def first_last_cancel(ls):
    return (ls[0] == 0 and ls[-1] != 0) or (ls[-1] == 0 and ls[0] != 0)

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

# Return all first/last-term-canceling partitions for correction n
def smart_partition(n):
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

def combine_first_last(ls):
    if first_last_cancel(ls):
        raise InvalidOperation()
    if ls[0] == 0:
        return DeltaTerm(-1, ls[:-1])
    return DeltaTerm(1, [ls[0] + ls[-1]] + ls[1:-1])

class DeltaTerm:
    def __init__(self, coeff, ls):
        self.coeff = coeff
        self.ls = ls

    def __eq__(self, other) -> bool:
        return cyclic_list_equal(other.ls, self.ls) and other.coeff == self.coeff

    def __add__(self, other):
        if not cyclic_list_equal(other.ls, self.ls):
            raise InvalidOperation()
        return DeltaTerm(self.coeff + other.coeff)

    def __str__(self) -> str:
        return "{:d} * {}".format(self.coeff, self.ls)

    def __repr__(self) -> str:
        return str(self)

def resolve_signs(dterm):
    return DeltaTerm(dterm.coeff * ((-1)**sum([1 for i in range(len(dterm.ls)) if dterm.ls[i] == 0])), dterm.ls)

def lists_equal(l1, l2):
    if len(l1) != len(l2):
        return False
    for i in range(len(l1)):
        if l1[i] != l2[i]:
            return False
    return True

class MultiSet:
    def __init__(self) -> None:
        self.items = []

    def add(self, new_item, count):
        need_to_append = True
        for i in range(len(self.items)):
            current_item, current_count = self.items[i]
            if current_item == new_item:
                self.items[i] = (current_item, current_count + count)
                need_to_append = False
                break
        if need_to_append:
            self.items.append(tuple((new_item, count)))
        self.clear_zero_count_items()

    def clear_zero_count_items(self):
        self.items = [i for i in self.items if i[1] != 0]

    def __eq__(self, other):
        if len(self.items) != len(other.items):
            return False
        for item in self.items:
            val, count = item
            match_found = False
            for i in other.items:
                if i[0] == val and i[1] == count:
                    match_found = True
                    break
            if not match_found:
                return False
        return True

    def __str__(self):
        return " ".join(["{}^{:d}".format(i[0], i[1]) if i[1] != 1 else "{}".format(i[0]) for i in self.items])

    def __repr__(self):
        return str(self)

class SigmaTerm:
    def __init__(self,exps):
        self.exps = exps

    def __eq__(self, other: object) -> bool:
        return lists_equal(self.exps, other.exps)
    
    def __str__(self):
        return "Sigma({})".format(self.exps)

    def __repr__(self) -> str:
        return str(self)

    def ToLatex(self):
        return "\\Sigma_{{}}".format(self.exps)

class PerturbativeTerm:
    def __init__(self, v_exp, sigmas, coeff) -> None:
        self.v_exp = v_exp
        self.sigmas = sigmas
        self.coeff = coeff

    def __eq__(self, other):
        return self.v_exp == other.v_exp and self.coeff == other.coeff and self.sigmas == other.sigmas

    def __str__(self) -> str:
        return "{:d} * (V_00 ^ {:d}) ".format(self.coeff, self.v_exp) + str(self.sigmas) if self.v_exp > 0 else "{:d} * ".format(self.coeff) + str(self.sigmas)
    
    def __repr__(self) -> str:
        return str(self)

def DeltaToPerturb(dterm):
    if dterm.ls[0] != 0:
        #rotate to start at a 0
        i = dterm.ls.index(0)
        dterm.ls = dterm.ls[i:] + dterm.ls[:i]
    
    pterm = PerturbativeTerm(0, MultiSet(), dterm.coeff)
    
    # get info for v_exp
    zero_run_lens = [1]
    on_run = True
    for i in range(1, len(dterm.ls)):
        if dterm.ls[i] == 0:
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
    for i in range(1, len(dterm.ls)):
        if dterm.ls[i] != 0:
            if in_group:
                number_groups[-1].append(dterm.ls[i])
            else:
                number_groups.append([dterm.ls[i]])
                in_group = True
        else:
            in_group = False
    for group in number_groups:
        pterm.sigmas.add(SigmaTerm(group), 1)
    
    return pterm

def print_pterms_by_v00(pterms):
    v_max = max(pterms, key=lambda t: t.v_exp).v_exp
    for i in range(v_max + 1):
        for t in [t for t in pterms if t.v_exp == i]:
            print(t)
    
n = 6
possible_terms = smart_partition(n)
combined_terms = [combine_first_last(ls) for ls in possible_terms] # combine the first and last indicies of each term, creating list of DeltaTerm objects
resolved_terms = [resolve_signs(ls) for ls in combined_terms] # resolve the signs of the DeltaTerm objects (odd # of 0s add a - sign)

# Consolidate the terms together (add like terms)
consolidated_terms = [resolved_terms[0]]
for term in resolved_terms[1:]:
    need_to_append = True
    for inserted_term in consolidated_terms:
        if cyclic_list_equal(term.ls, inserted_term.ls):
            inserted_term.coeff += term.coeff
            need_to_append = False
            break
    if need_to_append:
        consolidated_terms.append(term)

consolidated_terms = [t for t in consolidated_terms if t.coeff != 0] # remove terms w/ coeff 0
pterms = [DeltaToPerturb(term) for term in consolidated_terms] # convert from DeltaTerms to PerturbativeTerms

# Consolidate pterms together (add like terms)
consolidated_pterms = [pterms[0]]
for term in pterms[1:]:
    need_to_append = True
    for inserted_term in consolidated_pterms:
        if term.sigmas == inserted_term.sigmas:
            inserted_term.coeff += term.coeff
            need_to_append = False
            break
    if need_to_append:
        consolidated_pterms.append(term)

consolidated_pterms = [t for t in consolidated_pterms if t.coeff != 0] # remove terms w/ coeff 0
#print_pterms_by_v00(consolidated_pterms)
print(len(consolidated_pterms))