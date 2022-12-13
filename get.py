import requests
import json


try:
    with open('creds.txt', mode='r', encoding='UTF-8') as file:
        creds = json.load(file)
except FileNotFoundError:
    print('please create a file named creds.txt with the following format:')
    print('{"user": "your_username_or_email", "passcode": "your_passcode"}')
    exit(1)
except json.decoder.JSONDecodeError:
    creds = {}

if not 'user' in creds or not 'passcode' in creds:
    print('please make sure creds.txt is in the following format:')
    print('{"user": "your_username_or_email", "passcode": "your_passcode"}')
    exit(1)

try:
    cookies = requests.post('https://www.space-track.org/ajaxauth/login', data={'identity': creds['user'], 'password': creds['passcode']}).cookies
except requests.exceptions.ConnectionError:
    print('please check your internet connection')
    exit(1)
except requests.exceptions.HTTPError:
    print('please check your credentials')
    exit(1)

js: dict | list | str | int | float | bool | None = json.loads(requests.get('https://www.space-track.org/basicspacedata/query/class/gp/OBJECT_TYPE/Debris/MEAN_MOTION/12.2--15.6/', cookies=cookies).text)
if 'error' in js:
    print('please ensure your credentials are correct')
    exit(2)

with open('data/leo_debris.json', mode='w', encoding='UTF-8') as file:
    json.dump(js, file, indent=4)
