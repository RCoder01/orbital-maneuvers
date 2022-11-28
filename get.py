import requests
import json

with open('creds.txt', mode='r', encoding='UTF-8') as file:
    cookies = requests.post('https://www.space-track.org/ajaxauth/login', data={'identity': file.readline().replace('\n', ''), 'password': file.readline().replace('\n', '')}).cookies

js: dict | list | str | int | float | bool | None = json.loads(requests.get('https://www.space-track.org/basicspacedata/query/class/gp/OBJECT_TYPE/Debris/MEAN_MOTION/12.2--15.6/', cookies=cookies).text)
with open('data/leo_debris.json', mode='w', encoding='UTF-8') as file:
    json.dump(js, file, indent=4)
