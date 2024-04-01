from NDPTE import EnergyCorrection
from TexSoup import TexSoup
from datetime import datetime
import pickle
import numpy as np
import re

def string_count(string, substrings):
    count = 0
    for i in range(len(string)):
        for s in substrings:
            if i + len(s) <= len(string) and string[i:i+len(s)] == s:
                count += 1
    return count

def sigma_factor_to_latex(sf, exp=1) -> str:
        if exp == 0:
            return ""
        elif exp == 1:
            return r"\Sigma_{{{ind}}}".format(ind=",".join([str(x) for x in sf.indices]))
        else:
            return r"\Sigma_{{{ind}}}^{{{exp}}}".format(ind=",".join([str(x) for x in sf.indices]), exp=exp)
        
def perturbative_term_to_latex(pterm, coeff, with_V00=True) -> str:
        if coeff == 0:
            return ""
        if abs(coeff) > 1:
            coeff_str = str(int(coeff))
            if coeff > 0:
                coeff_str = "+" + coeff_str
        else:
            coeff_str = "+" if coeff==1 else "-"
        sigmas = "".join([sigma_factor_to_latex(sigma, exp) for sigma, exp in pterm.sigmas.items()])
        if with_V00 and pterm.v_exp != 0:
            raise NotImplementedError
        else:
            return coeff_str + sigmas

def line_length(line):
    last_line = line[line.rfind(r"\\"):] if r"\\" in line else line
    line1 = re.sub("\^\{[\d]+\}", "", last_line) 
    line1 = re.sub("_\{[\d,]+\}", "", line1)
    line2 = re.sub("\\Sigma", "S", line1)
    line3 = line2.replace("{", "").replace("}", "").replace("\\", "").replace(" ", "").replace("\n\t\t&quadquad", "")
    s_count = line3.count("S")
    pm_count = string_count(line3, ["+", "-"])
    num_count = string_count(line3, ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"])
    underscore_num_count = 0
    underscore_comma_count = 0
    for s in re.findall("_\{[\d]+\}", last_line):
        underscore_num_count += string_count(s, ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"])
        underscore_comma_count += s.count(",")
    foo = s_count + 1.5*pm_count + 0.8*num_count + 0.6*underscore_num_count + 0.15*underscore_comma_count
    return foo

def correction_to_latex(ec):
        res = "\t\\trdelta{" + str(ec.n) + "} &= "
        for v_exp, terms in sorted(ec.TermsByV00.items(), key=lambda x:x[0]):
            gcd = np.gcd.reduce([t[1] for t in terms])
            neg = len([t for t in terms if t[1] < 0]) > len(terms) / 2
            coeff_str = str(int(gcd)) if gcd != 1 else ""
            coeff_str = ("-" if neg else "+") + coeff_str
            sigma_terms = ""
            sigma_and_subscript_count = 0
            newline = False
            first_line = True
            for term, coeff in terms:    
                sigma_and_subscript_count += sum([1 + len(sigma.indices)for sigma in term.sigmas.elements()])
                    
                sigma_terms += perturbative_term_to_latex(term, coeff/gcd if not neg else -coeff/gcd, with_V00=False)
                line_len = line_length(sigma_terms) if not first_line else line_length(sigma_terms) + 4
                newline = line_len > 20
                if newline:
                    first_line = False
                    sigma_terms += "({:.3f})".format(line_len) +  r" \\ " + "\n\t\t&\\quad\\quad "
                    sigma_and_subscript_count = 0
            if newline:
                sigma_terms = sigma_terms[:-18]
            if sigma_terms[0] == "+":
                sigma_terms = sigma_terms[1:]
            if v_exp > 0:
                v_str = r"V_{00}"
                if v_exp > 1:
                    v_str += "^{" + (str(v_exp) if v_exp > 1 else "") + "}"
                res += "\t\t&\\quad\\; " + coeff_str + v_str
            else:
                res += (coeff_str if not coeff_str in ["+", "-"] else "")
            if len(terms) > 1 and v_exp > 0:
                res += r"(" + sigma_terms + ")"
            else:
                res += sigma_terms
            res += r" \\" + "\n"
        return res

with open("original_tex.tex", 'r') as tex_file:
    tex = TexSoup(tex_file)
align = tex.find("align*")

with open("data/corrections.pkl", 'rb') as data_file:
    data = pickle.load(data_file)

for correction_data in data.items():
    if correction_data[0] != 12:
        continue
    start_time = datetime.now()
    print(correction_data[0])
    print("Start Time: ", start_time)

    correction = EnergyCorrection.from_tuple(correction_data)
    correction_str = correction_to_latex(correction)
    correction_tex = TexSoup(correction_str)

    align.append(correction_tex)

    with open("output/wtf.tex", "w") as f:
        f.write(str(tex))
    
    end_time = datetime.now()
    print("End Time: ", end_time)
    diff = end_time - start_time
    print("Elapsed Time: {} days {} hours {} minutes {} seconds {} milliseconds {} microseconds".format(diff.days, diff.seconds // 3600, (diff.seconds // 60) % 60, diff.seconds % 60, diff.microseconds // 1000, diff.microseconds % 1000))
    print("Characters: {}".format(len(correction_str)))
    print("Lines: {}".format(correction_str.count("\n")))
    print()