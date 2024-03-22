from NDPTE import EnergyCorrection
from TexSoup import TexSoup

with open("original_tex.tex", 'r') as tex_file:
    tex = TexSoup(tex_file)
dgroup = tex.find("dgroup*")

for i in range(4, 10):
    print(i)
    correction = EnergyCorrection(i)
    correction.calc()
    correction.sort_by_v00()

    correction_str = correction.to_latex()
    correction_tex = TexSoup(correction_str)

    dgroup.append(correction_tex)

with open("output/energies.tex", "w") as f:
    f.write(str(tex))