#!/usr/bin/env python
"""
Script to pull data, and process through velocity results and time series time series plots
"""

import requests
from pathlib import Path
 
from earthscope_sdk.auth.device_code_flow import DeviceCodeFlowSimple
from earthscope_sdk.auth.auth_flow import NoTokensError

# choose where you want the token saved - the default file name is sso_tokens.json
# if you want to keep the default name, set the path to a directory. Include a file name to rename. 
token_path = "/Users/an77/tokens/"

url = "https://data-idm.unavco.org/path/to/data_file"

url = "https://data.unavco.org/archive/gnss/products/position/BYAH/BYAH.cwu.igs14.csv"
# example: "https://data-idm.unavco.org/archive/gnss/rinex/obs/2022/298/ar272980.22d.Z"

# instantiate the device code flow subclass
device_flow = DeviceCodeFlowSimple(Path(token_path))
try:
    # get access token from local path
    device_flow.get_access_token_refresh_if_necessary()
except NoTokensError:
    # if no token was found locally, do the device code flow
    device_flow.do_flow()
token = device_flow.access_token
 
# request a file and provide the token in the Authorization header
file_name = Path(url).name
directory_to_save_file = Path.cwd() # where you want to save the downloaded file 

r = requests.get(url, headers={"authorization": f"Bearer {token}"})
if r.status_code == requests.codes.ok:
    # save the file
    with open(Path(directory_to_save_file / file_name), 'wb') as f:
        for data in r:
            f.write(data)
else:
    #problem occured
    print(f"failure: {r.status_code}, {r.reason}")