# -*- coding: utf-8 -*-
"""
Sending a set of commands to a Bruker 125HR and saving the returned opus file

Author: Matthias Buschmann, IUP Bremen
Date: 2023/06/01
"""

import os, sys, requests

def read_commands(fname):
    with open(fname, 'r', encoding='utf8') as f:
        ll = f.readlines()
    return [l.strip() for l in ll if len(l.strip())>1 and not l.strip().startswith('#')]
    
def save_opus(fname, o):
    with open(fname, 'wb') as f:
        f.write(o)

def get_status(stat_htm):
    stat = requests.get()
    i1 = stat.text.rfind('ID=MSTCO')
    i2 = stat.text.find('<', i1)
    status = stat.text[i1 + 9:i2]
    return status

def run_measurement(url_ftir, meas_htm, data_htm, stat_htm):
    requests.get(meas_htm)
    status = 'SCN'
    data = None
    while status != 'IDL':
        # repeat requests until IDL
        status = get_status()
        if status=='IDL':
            # find download link
            data = requests.get(data_htm)
            i1 = data.text.find('A HREF=')
            i2 = data.text.find('">', i1)
            # download data from ifs
            data = requests.get('/'.join((url_ftir,data.text[i1+9:i2])))
            datac = data.content
            break
    return datac    


if __name__ == '__main__':
    if len(sys.argv!=4):
        print('Usage:\n\t python record_dip_per_apt.py path/to/commandlist path/to/save/to/ filename_prefix \n\ne.g.: python record_dip_per_apt.py apt_test_ingaas_br data/br_ingaas_commands data/ br_ingaas_apt_test_')
    cmd_file = sys.argv[1]
    save_folder = sys.argv[2]
    fname_prefix = sys.argv[3]
    cmds = read_commands(cmd_file)
    ip = cmds[0]
    url_ftir = 'http://'+ip
    stat_htm = '/'.join((url_ftir,'stat.htm'))
    data_htm = '/'.join((url_ftir, 'datafile.htm'))
    for i,cmd in enumerate(cmds[1:]):
        meas_htm = '/'.join((url_ftir, cmd))
        print('Running command: ', meas_htm)
        o = run_measurement(url_ftir, meas_htm, data_htm, stat_htm)
        fname = os.path.join(savefolder, prefix+f'{i:02d}')
        save_opus(fname, o)
        print('Saved fname')