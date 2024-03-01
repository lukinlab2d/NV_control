"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""
# Python program to read
# json file
  
  
import json
  
# Opening JSON file
# C:\Users\lukin2dmaterials\data\2023-04-07\#020_ODMR_CW_19-45-29
# filename = 'C:/Users/lukin2dmaterials/data/2023-04-12/#055_T2E_17-52-52/snapshot.json'
filename = 'C:/Users/lukin2dmaterials/data/2023-06-08/#043_CalibReadoutPhotonStat_18-15-15/snapshot.json'
f = open(filename)
  
# returns JSON object as 
# a dictionary
data = json.load(f)
  
# Iterating through the json list
for k, v in data.items():
    if k != 'loop' and k!= 'arrays':
        print (k, v)

# for i in data['loop']:
#     print(i)
    
# print('---------')
# for ii in data['loop']['sweep_values']:
#     print(ii)
#     print(data['loop']['sweep_values'][ii])

# print('---------')
# print( data['loop']['delay'])

# print('---------')
# print(data['loop']['use_threads'])
  
# # Closing file
# f.close()