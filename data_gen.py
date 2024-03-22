from NDPTE import EnergyCorrection
import pickle

for i in range(13,21):
    correction = EnergyCorrection(i)
    correction.calc()
    correction.sort_by_v00()
    data = correction.to_tuple()
    print("Terms in correction {}: {}".format(i, len(correction.p_terms.items())))
    with open("data/Corrections13-20.pkl", 'ab') as file:
        pickle.dump(data, file)