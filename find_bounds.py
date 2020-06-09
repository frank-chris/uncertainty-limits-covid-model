import os
import pandas as pd
import numpy as np

def statename(statecode):
    '''
    Returns state name from state code   
    '''
    if statecode == "AP":
        return "Andhra Pradesh"
    elif statecode == "AN":
        return "Andaman and Nicobar Islands"
    elif statecode == "AR":
        return "Arunachal Pradesh"
    elif statecode == "AS":
        return "Assam"
    elif statecode == "BR":
        return "Bihar"
    elif statecode == "CH":
        return "Chandigarh"
    elif statecode == "CT":
        return "Chhattisgarh"
    elif statecode == "DD":
        return "Dadra and Nagar Haveli and Daman and Diu"
    elif statecode == "DL":
        return "Delhi"
    elif statecode == "GA":
        return "Goa"
    elif statecode == "GJ":
        return "Gujarat"
    elif statecode == "HR":
        return "Haryana"
    elif statecode == "HP":
        return "Himachal Pradesh"
    elif statecode == "JH":
        return "Jharkhand"
    elif statecode == "JK":
        return "Jammu and Kashmir"
    elif statecode == "KA":
        return "Karnataka"
    elif statecode == "KL":
        return "Kerala"
    elif statecode == "LA":
        return "Ladakh"
    elif statecode == "LD":
        return "Lakshadweep"
    elif statecode == "MP":
        return "Madhya Pradesh"
    elif statecode == "MH":
        return "Maharashtra"
    elif statecode == "MN":
        return "Manipur"
    elif statecode == "ML":
        return "Meghalaya"
    elif statecode == "MZ":
        return "Mizoram"
    elif statecode == "NL":
        return "Nagaland"
    elif statecode == "OR":
        return "Odisha"
    elif statecode == "PB":
        return "Punjab"
    elif statecode == "PY":
        return "Puducherry"
    elif statecode == "RJ":
        return "Rajasthan"
    elif statecode == "SK":
        return "Sikkim"
    elif statecode == "TN":
        return "Tamil Nadu"
    elif statecode == "TG":
        return "Telengana"
    elif statecode == "TR":
        return "Tripura"
    elif statecode == "UP":
        return "Uttar Pradesh"
    elif statecode == "UT":
        return "Uttarakhand"
    elif statecode == "WB":
        return "West Bengal"
    else:
        return statecode

def cumulative(df):
    for column in df.columns:
        df[column] = df[column].cumsum()
    return df


d = '.'
folders = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]
folders = [x.replace('.\\', '') for x in folders]
print('\n'+ str(len(folders)) + ' folders found.\n')

folders_with_deaths = [x for x in folders if 'CovidDeaths.data' in os.listdir(x+'/')]

print('\n'+ str(len(folders_with_deaths)) + ' folders with CovidDeaths.data found.\n')

df = pd.read_csv(folders[0]+'/CovidPopulation.data', delimiter=" ", header=1)
del df['Day']
states = list(df.columns)
indices = list(df.index)

min_df = pd.read_csv(folders[0]+'/CovidPopulation.data', delimiter=" ", header=1)
max_df = pd.read_csv(folders[0]+'/CovidPopulation.data', delimiter=" ", header=1)

for i in range(1,len(folders)):
    df = pd.read_csv(folders[i]+'/CovidPopulation.data', delimiter=" ", header=1)
    for state in states:
        for index in indices:
            if df.loc[index, state] > max_df.loc[index, state]:
                max_df.loc[index, state] = df.loc[index, state]
            if df.loc[index, state] < min_df.loc[index, state]:
                min_df.loc[index, state] = df.loc[index, state]
            
hightotalData = {}
lowtotalData = {}

states.remove('Total')

states.append('DD')
states.append('ML')
states.append('MZ')
states.append('NL')

highstatesData = {}
lowstatesData = {}

for index in indices:
    hightotalData[str(index)] = max_df.loc[index,'Total']
    lowtotalData[str(index)] = min_df.loc[index,'Total']

min_df['DD'] = 0
min_df['ML'] = 0
min_df['MZ'] = 0
min_df['NL'] = 0
max_df['DD'] = 0
max_df['ML'] = 0
max_df['MZ'] = 0
max_df['NL'] = 0





for state in states:
    highstatesData[statename(state)] = {}
    lowstatesData[statename(state)] = {}
    for index in indices:
        highstatesData[statename(state)][str(index)] = max_df.loc[index, state]
        lowstatesData[statename(state)][str(index)] = min_df.loc[index, state]



df = pd.read_csv(folders[0]+'/CovidRecovered.data', delimiter=" ", header=1)
del df['Day']
states = list(df.columns)
indices = list(df.index)

min_df = pd.read_csv(folders[0]+'/CovidRecovered.data', delimiter=" ", header=1)
min_df = cumulative(min_df)
max_df = pd.read_csv(folders[0]+'/CovidRecovered.data', delimiter=" ", header=1)
max_df = cumulative(max_df)


for i in range(1,len(folders)):
    df = pd.read_csv(folders[i]+'/CovidRecovered.data', delimiter=" ", header=1)
    df = cumulative(df)
    for state in states:
        for index in indices:
            if df.loc[index, state] > max_df.loc[index, state]:
                max_df.loc[index, state] = df.loc[index, state]
            if df.loc[index, state] < min_df.loc[index, state]:
                min_df.loc[index, state] = df.loc[index, state]


for index in indices:
    hightotalData['Recovered'+str(index)] = max_df.loc[index,'Total']
    lowtotalData['Recovered'+str(index)] = min_df.loc[index,'Total']


states.remove('Total')
states.append('DD')
states.append('ML')
states.append('MZ')
states.append('NL')


min_df['DD'] = 0
min_df['ML'] = 0
min_df['MZ'] = 0
min_df['NL'] = 0
max_df['DD'] = 0
max_df['ML'] = 0
max_df['MZ'] = 0
max_df['NL'] = 0



for state in states:
    for index in indices:
        highstatesData[statename(state)]['Recovered'+str(index)] = max_df.loc[index, state]
        lowstatesData[statename(state)]['Recovered'+str(index)] = min_df.loc[index, state]

# deceased

df = pd.read_csv(folders_with_deaths[0]+'/CovidDeaths.data', delimiter=" ", header=1)
del df['Day']
states = list(df.columns)
indices = list(df.index)

min_df = pd.read_csv(folders_with_deaths[0]+'/CovidDeaths.data', delimiter=" ", header=1)
min_df = cumulative(min_df)
max_df = pd.read_csv(folders_with_deaths[0]+'/CovidDeaths.data', delimiter=" ", header=1)
max_df = cumulative(max_df)


for i in range(1,len(folders_with_deaths)):
    df = pd.read_csv(folders_with_deaths[i]+'/CovidDeaths.data', delimiter=" ", header=1)
    df = cumulative(df)
    for state in states:
        for index in indices:
            if df.loc[index, state] > max_df.loc[index, state]:
                max_df.loc[index, state] = df.loc[index, state]
            if df.loc[index, state] < min_df.loc[index, state]:
                min_df.loc[index, state] = df.loc[index, state]


for index in indices:
    hightotalData['Deceased'+str(index)] = max_df.loc[index,'Total']
    lowtotalData['Deceased'+str(index)] = min_df.loc[index,'Total']


states.remove('Total')
states.append('DD')
states.append('ML')
states.append('MZ')
states.append('NL')


min_df['DD'] = 0
min_df['ML'] = 0
min_df['MZ'] = 0
min_df['NL'] = 0
max_df['DD'] = 0
max_df['ML'] = 0
max_df['MZ'] = 0
max_df['NL'] = 0



for state in states:
    for index in indices:
        highstatesData[statename(state)]['Deceased'+str(index)] = max_df.loc[index, state]
        lowstatesData[statename(state)]['Deceased'+str(index)] = min_df.loc[index, state]


with open("uncertainty.js", 'w') as file:
    file.write("var highstatesData = " + str(highstatesData) + ";"+"var lowstatesData = " + str(lowstatesData) + ";"+"var hightotalData = " + str(hightotalData) + ";" + "var lowtotalData = " + str(lowtotalData) + ";" )
print("\nData written into uncertainty.js")

