# uncertainty-limits-covid-model

## Description

### 1. find_bounds.py

Run this script to write the extreme values(highs and lows) into **uncertainty.js.** find_bounds.py 
uses the data from all the folders present in the repository(there are 17 such folders now) and writes the lowest and highest values into uncertainty.js.
Note that it creates uncertainty for deaths from only those folders that have a CovidDeaths.data file(5 such folders in the total 17 folders).

### 2. folders

The folders contain the data which is used for computing the uncertainty.

### 3. chart.js, data.js, index.html, plotactivecases.py

These are the files created when I was asked to plot a few runs(just before the first time uncertainty was plotted) to help choose among the runs.
Open the html file to see the plots.

### 4. uncertainty.js

The python script named find_bounds.py writes the computed uncertainty into uncertainty.js. 
This file should replace the uncertainty.js file in the webpage's repo(in folder named data/uncertainty), for changes to take place in the webpage.

## Usage

1. Keep only those folders that should be considered for computing the uncertainty. 

2. Run find_bounds.py

3. Copy the uncertainty.js file and replace the previous one in the webpage repo with this newly created one.
