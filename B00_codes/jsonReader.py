# Python program to read
# json file
  
  
import json
  
# Opening JSON file
filename = 'C:/Users/lukin2dmaterials/data/2022-12-13/#008_T2E_02-26-33/snapshot.json'
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