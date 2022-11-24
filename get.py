import requests
import json

with open('creds.txt', mode='r', encoding='UTF-8') as f:
    cookies = requests.post('https://www.space-track.org/ajaxauth/login', data={'identity': f.readline().replace('\n', ''), 'password': f.readline().replace('\n', '')}).cookies

js: dict | list | str | int | float | bool | None = json.loads(requests.get('https://www.space-track.org/basicspacedata/query/class/gp/OBJECT_TYPE/Debris/MEAN_MOTION/12.2--15.6/', cookies=cookies).text)
