# For calculating energy corrections in Nondegenerate perturbation theory

import numpy as np
from collections import defaultdict
from typing import Callable, Generator
import pickle

# rotate list ls to the right by n (negative to go left)
def rotate_list(ls: list, n: int) -> list:
    return ls[-n:] + ls[:-n]

# Determine if two lists are rotated versions of each other
def list_cyclic_equal(l1: list, l2: list) -> list:
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
def list_set_equal(l1: list, l2: list) -> list:
    return np.all(np.array([x in l2 for x in l1])) and np.all(np.array(x in l1 for x in l2))

class PartitionGenerator:
    # Store previously calculated partitions for reuse
    partitions = dict()
    
    # Generate all lists of m non-negative integers that add to n
    def partition(self, n: int, m: int) -> list:
        if (n,m) in PartitionGenerator.partitions.keys():
            return PartitionGenerator.partitions[(n,m)]
        res = []
        if m == 0:
            return [[0 for _ in range(n)]]
        if n == 1:
            return [[m]]
        for i in range(m+1):
            for ls in self.partition(n-1, m-i):
                res += [[i] + ls]
        PartitionGenerator.partitions[(n,m)] = res
        return res
    
    # Generates partitions of
    def noncanceling_partitions(self, n: int) -> Generator[list, None, None]:
        # zero on each side
        for ls in self.partition(n-1, n-1):
            yield [0] + ls + [0]
        # nonzero on each side
        for i in range(1, n):
            for j in range(1, n-i):
                for ls in self.partition(n-1, n-1-i-j):
                    yield [i] + ls + [j]

class MultiSet:
    def __init__(self, sort_key: Callable = lambda x:x) -> None:
        self.sort_key = sort_key
        self._dict = defaultdict(int)

    def add(self, new_item, count: int, auto_clear: bool = True) -> None:
        self._dict[new_item] += count
        if auto_clear:
            self.clear_zero_count_items()

    def add_many(self, items, counts: list[int]) -> None:
        for item, count in zip(items, counts):
            self.add(item, count, auto_clear=False)
        self.clear_zero_count_items()

    def clear_zero_count_items(self) -> None:
        new_items = defaultdict(int)
        for k,v in self._dict.items():
            if v != 0:
                new_items[k] = v
        self._dict = new_items

    def items(self):
        return self._dict.items()
    
    def elements(self):
        return self._dict.keys()
    
    def counts(self):
        return self._dict.values()

    def __eq__(self, other):
        return self._dict == other._dict
    
    def __hash__(self) -> int:
        return hash(tuple(sorted(self.items(), key=lambda x: self.sort_key(x[0]))))

    def __str__(self):
        return " ".join(["{}^{:d}".format(i[0], i[1]) if i[1] != 1 else "{}".format(i[0]) for i in self._dict.items()])

    def __repr__(self):
        return str(self)

class SigmaFactor:
    def __init__(self, indices: list[int]):
        self.indices = indices

    def __eq__(self, other: object) -> bool:
        return self.indices == other.indices
    
    def __hash__(self):
        return hash(tuple(self.indices))

    def __str__(self) -> str:
        return "Sigma({})".format(self.indices)

    def __repr__(self) -> str:
        return str(self)
    
    def to_latex(self, exp=1) -> str:
        if exp == 0:
            return ""
        elif exp == 1:
            return r"\Sigma_{{{ind}}}".format(ind=",".join([str(x) for x in self.indices]))
        else:
            return r"\Sigma_{{{ind}}}^{{{exp}}}".format(ind=",".join([str(x) for x in self.indices]), exp=exp)

class PerturbativeTerm:
    def __init__(self, v_exp: int, sigmas) -> None:
        self.v_exp = v_exp
        self.sigmas = sigmas # Multiset of SigmaFactors

    @staticmethod
    def FromPartition(partition: list[int]):
        # Combine first/last
        partition = [partition[0] + partition[-1]] + partition[1:-1]

        #rotate to start at a boundary 0
        for i in range(len(partition)):
            if partition[i] != 0 and partition[(i+1)%len(partition)] == 0:
                shift = i+1    
        partition = rotate_list(partition, -shift)
        
        pterm = PerturbativeTerm(0, MultiSet(sort_key=lambda x: str(x.indices)))
        
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
        pterm.sigmas.add_many(*zip(*[(SigmaFactor(group), 1) for group in number_groups]))
        
        return pterm

    def __eq__(self, other) -> bool:
        return self.v_exp == other.v_exp and self.sigmas == other.sigmas
    
    def __hash__(self):
        return hash((self.v_exp, self.sigmas))

    def __str__(self) -> str:
        return "(V_00 ^ {:d}) ".format(self.v_exp) + str(self.sigmas) if self.v_exp > 0 else ""
    
    def string_with_coeff(self, coeff) -> str:
        return "{:d} * (V_00 ^ {:d}) ".format(coeff, self.v_exp) + str(self.sigmas) if self.v_exp > 0 else "{:d} * ".format(coeff) + str(self.sigmas)

    def to_latex(self, coeff, with_V00=True) -> str:
        if coeff == 0:
            return ""
        if abs(coeff) > 1:
            coeff_str = str(int(coeff))
            if coeff > 0:
                coeff_str = "+" + coeff_str
        else:
            coeff_str = "+" if coeff==1 else "-"
        sigmas = "".join([sigma.to_latex(exp) for sigma, exp in self.sigmas.items()])
        if with_V00 and self.v_exp != 0:
            raise NotImplementedError
        else:
            return coeff_str + sigmas

    def __repr__(self) -> str:
        return str(self)

class EnergyCorrection:
    def __init__(self, n):
        self.n = n
        self.TermsByV00 = dict()
        self.p_terms = MultiSet()

    def calc(self):
        # Generate perms and combine first/last 
        p_terms = MultiSet() # Multiset of SigmaFactors
        partGen = PartitionGenerator()
        for part in partGen.noncanceling_partitions(self.n):
            p_terms.add(PerturbativeTerm.FromPartition(part), (-1)**(part.count(0)), auto_clear=False)
        p_terms.clear_zero_count_items()
        self.p_terms = p_terms

    def sort_by_v00(self):
        #Sort by V00
        self.v_max = max(list(self.p_terms.elements()), key=lambda t: t.v_exp).v_exp
        for i in range(self.v_max + 1):
            self.TermsByV00[i] = [t for t in self.p_terms.items() if t[0].v_exp == i]

    def print_pterms_by_v00(self):
        for v_exp, terms in self.TermsByV00.items():
            for t in terms:
                print(t[0].string_with_coeff(t[1]))
        
    def to_tuple(self):
        return (self.n, [(v_exp, [([(sigma.indices, exp) for sigma,exp in sigmas.sigmas.items()], coeff) for sigmas,coeff in pterm]) for v_exp, pterm in self.TermsByV00.items()])
    
    def to_latex(self):
        res = "\t\\begin{dmath*}\n\t\t\\trdelta{" + str(self.n) + "} = "
        for v_exp, terms in sorted(self.TermsByV00.items(), key=lambda x:x[0]):
            gcd = np.gcd.reduce([t[1] for t in terms])
            neg = len([t for t in terms if t[1] < 0]) > len(terms) / 2
            coeff_str = str(int(gcd)) if gcd != 1 else ""
            coeff_str = ("-" if neg else "+") + coeff_str
            sigma_terms = ""
            for term, coeff in terms:
                sigma_terms += term.to_latex(coeff/gcd if not neg else -coeff/gcd, with_V00=False)
            if sigma_terms[0] == "+":
                sigma_terms = sigma_terms[1:]
            if v_exp > 0:
                v_str = r"V_{00}"
                if v_exp > 1:
                    v_str += "^{" + (str(v_exp) if v_exp > 1 else "") + "}"
                res += "\t\t\\quad " + coeff_str + v_str
            else:
                res += (coeff_str if not coeff_str in ["+", "-"] else "")
            if len(terms) > 1 and v_exp > 0:
                res += r"(" + sigma_terms + ")"
            else:
                res += sigma_terms
            res += r" \\" + "\n"
        res += "\t\\end{dmath*}\n"
        return res

    @staticmethod
    def from_tuple(tup):
        n, termsByVexp = tup
        correction = EnergyCorrection(n)
        for v_exp, terms in termsByVexp:
            for term,coeff in terms:
                sigmas = MultiSet(sort_key=lambda x: str(x.indices))
                for ind,exp in term:
                    sigmas.add(SigmaFactor(ind), exp)
                correction.p_terms.add(PerturbativeTerm(v_exp, sigmas), coeff)
        correction.sort_by_v00()
        return correction