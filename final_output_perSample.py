#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 15:32:20 2020

@author: User
"""


import pandas as pd

pdata = pd.read_csv("all_info_filtered_AD-DP.tab", delimiter="\t")

df = pdata

list_of_names = df['SAMPLE'].values.tolist()
list_of_calls = df['REC_CALL'].values.tolist()


d = {}

for name, pro in zip(list_of_names, list_of_calls):
    if name not in d:
        d[name] = pro
    elif type(d[name]) == list:
        d[name].append(pro)
    else:
        d[name] = [d[name], pro]


df_dict = pd.DataFrame(d.items(), columns = ["Sample", "GBA_variants"])

print(df_dict.head())

df_dict.to_csv(r'gba_pipeline_finalOutput.csv', index = False)
