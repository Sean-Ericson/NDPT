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

class SigmaTerm:
    def __init__(self, n, exps):
        self.n = n
        self.exps = exps

    def __eq__(self, other: object) -> bool:
        return self.n == other.n and lists_equal(self.exps, other.exps)
    
    def __str__(self):
        return "Sigma({:d}; {})".format(self.n, self.exps)

    def __repr__(self) -> str:
        return str(self)

    def ToLatex(self):
        return "\\Sigma_{{}}^{:d}".format(self.exps, self.n)

class PerturbativeTerm:
    def __init__(self, v_exp, sigmas, coeff) -> None:
        self.v_exp = v_exp
        self.sigmas = sigmas
        self.coeff = coeff

    def __str__(self) -> str:
        return "{:d} * (V_00 ** {:d}) ".format(self.coeff, self.v_exp) + " ".join([str(s) for s in self.sigmas]) if self.v_exp > 0 else "{:d} * ".format(self.coeff) + " ".join([str(s) for s in self.sigmas])
    
    def __repr__(self) -> str:
        return str(self)

class V0Term:
    def __init__(self, power, pterms) -> None:
        self.power = power
        self.pterms = pterms

    def __str__(self):
        overall_negative = len([x for x in self.pterms if x.coeff < 0]) >= math.ceil(len(self.pterms)/2)

    
def DeltaToPerturb(dterm):
    if dterm.ls[0] != 0:
        raise Exception("I thought (well, conjectured) I was always gonna have 0s at the start...")
    pterm = PerturbativeTerm(0, [], dterm.coeff)
    
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
        pterm.sigmas.append(SigmaTerm(len(group)-1, group))
    
    return pterm
    
n = 6
possible_terms = [ls for ls in partition(n+1,n-1) if not first_last_cancel(ls)] # only consider terms where the first/last indicies don't cancel
combined_terms = [combine_first_last(ls) for ls in possible_terms] # combine the first and last indicies of each term
resolved_terms = [resolve_signs(ls) for ls in combined_terms] # resolve negative signs

# Consolidate the terms together
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

pterms = [DeltaToPerturb(term) for term in consolidated_terms]

# Consolidate pterms together
consolidated_pterms = [pterms[0]]
for term in pterms[1:]:
    need_to_append = True
    for inserted_term in consolidated_pterms:
        if set_lists_equal(term.sigmas, inserted_term.sigmas):
            inserted_term.coeff += term.coeff
            need_to_append = False
            break
    if need_to_append:
        consolidated_pterms.append(term)


consolidated_pterms = [t for t in consolidated_pterms if t.coeff != 0]
#print(consolidated_pterms)
print([t for t in consolidated_pterms if t.v_exp == 1])
