import os
import pandas as pd

d = '.'
folders = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
folders = [x.replace('.\\', '') for x in folders]


df = pd.read_csv(folders[0]+'/CovidPopulation.data', delimiter=" ", header=1)
del df['Day']
states = list(df.columns)

data_dict = {}

for state in states:
    data_list = []
    for folder in folders:
        df = pd.read_csv(folder+'/CovidPopulation.data', delimiter=" ", header=1)
        for i in df.index:
            data_list.append([i, folder, df.loc[i, state]])
    data_dict[state] = data_list
    
with open("data.js", 'w') as file:
    file.write("var data = " + str(data_dict) + ";")

print("\nData written into data.js")

