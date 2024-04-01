from NDPTE import EnergyCorrection
import pickle
import matplotlib.pyplot as plt

with open("data/Corrections.pkl", 'rb') as file:
    data = pickle.load(file)


lens = {1: 1, 2: 1, 3: 2}
for dat in data.items():
    correction = EnergyCorrection.from_tuple(dat)
    count = correction.term_count()
    print("Terms in correction {:02d}: {:11,}".format(correction.n, count))
    lens[correction.n] = count

plt.plot(*zip(*list(lens.items())))
plt.yscale("log")
plt.title("Terms in nth Order NDPT Energy Correction")
plt.xlabel("Order")
plt.ylabel("Terms")
plt.show()