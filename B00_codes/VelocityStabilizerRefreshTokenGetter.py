"""
This file is part of B00 codes based on b26_toolkit. Questions are addressed to Hoang Le.
"""

import base64
import requests
import json

tkFile = 'C:/Users/lukin2dmaterials/data/tk.txt'; lines = []
with open(tkFile, 'r') as file:
    for line in file:
        if "\n" in line: line = line[0:-1]
        lines.append(line)

APP_KEY = lines[1]
APP_SECRET = lines[2]
ACCESS_CODE_GENERATED = 0

BASIC_AUTH = base64.b64encode(f'{APP_KEY}:{APP_SECRET}'.encode())

headers = {
    'Authorization': f"Basic {BASIC_AUTH.decode()}",
    'Content-Type': 'application/x-www-form-urlencoded',
}

data = f'code={ACCESS_CODE_GENERATED}&grant_type=authorization_code'

response = requests.post('https://api.dropboxapi.com/oauth2/token',
                         data=data,
                         auth=(APP_KEY, APP_SECRET))
print(json.dumps(json.loads(response.text), indent=2))