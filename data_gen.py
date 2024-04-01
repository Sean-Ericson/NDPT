from NDPTE import EnergyCorrection
import pickle
from datetime import datetime, timedelta

try:
    with open("data/Corrections.pkl", 'rb') as file:
        data = pickle.load(file)
    i = max(list(data.keys()))
except:
    i = 3
    data = {}

while True:
    i += 1
    print(i)
    start_time = datetime.now()
    print("Start Time: ", start_time)
    correction = EnergyCorrection(i)
    correction.calc()
    correction.sort_by_v00()
    d = correction.to_tuple()
    data[i] = d[1]
    
    with open("data/Corrections.pkl", 'wb') as file:
        pickle.dump(data, file)
    
    end_time = datetime.now()
    print("End Time: ", end_time)
    diff = end_time - start_time
    print("Elapsed Time: {} days {} hours {} minutes {} seconds {} milliseconds {} microseconds".format(diff.days, diff.seconds // 3600, (diff.seconds // 60) % 60, diff.seconds % 60, diff.microseconds // 1000, diff.microseconds % 1000))

    print("Terms in correction: {}".format(len(correction.p_terms.items())))
    print()