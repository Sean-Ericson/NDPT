from NDPTE import EnergyCorrection
import pickle

with open("data/Coded_Corrections.pkl", 'rb') as file:
    foo = pickle.load(file)
    while True:
        data = pickle.load(file)
        print(data)
        correction = EnergyCorrection.from_tuple(data)
        correction.print_pterms_by_v00()