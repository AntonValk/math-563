#!/usr/bin/python3

import os
import re

dir_path = os.path.dirname(os.path.realpath(__file__))

dirlist = None

def sorted_alphanumeric(data):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

for f in ['admm', 'chambollepock', 'primal_dr', 'primaldual_dr']:
        folder = dir_path + '/'+ f + '/' +  'Iterates'
        imgs = os.listdir(folder)
        dirlist = sorted_alphanumeric(os.listdir(folder))

        for i, name in enumerate(dirlist):
                words = name.split('_')
                new_name = '_'.join(words[:-1])+ '_iterate-' + str(i) + '.png'
                os.rename(folder + '/' + name, folder + '/' + new_name)

            
