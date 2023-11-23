import base64
import requests
import json

APP_KEY = 'jr02kdlisnsp67m'
APP_SECRET = 'ts6lhusxptmz4yk'
ACCESS_CODE_GENERATED = 'DMYteU69mToAAAAAAAAEEX921Yz0_likzbra8YkLeFA'

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