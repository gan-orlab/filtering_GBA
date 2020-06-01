#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 19:00:56 2020

@author: User
"""

import pandas as pd


def convert(list): 
      
    # Converting integer list to string list 
    s = [str(i) for i in list] 
      
    # Join list items using join() 
    res = "".join(s)
      
    return(res) 


def convert_gba_notation(variants):
    
    set_of_variants = set(variants)
    new_dict = {}
    
    for index, element in enumerate(set_of_variants):
        letters = []
        number = []
        garbage = []
        for i in element:
            if not i.isdecimal():
                letters.append(i)
            elif i.isdecimal():
                number.append(i)
            else: garbage.append(i)
        number = int(convert(number))
        if number > 39:
            new_number = number-39
        else:
            new_number = number-40
        
        new = [letters[0], str(new_number), letters[1]]
        new = convert(new)
        
        new_dict[element] = new
        
    
    df = pd.DataFrame(new_dict.items(), columns = ["PROTEIN_CHANGE", "NEW"])
    return df



df = pd.read_csv("likely_pathogenic_annotations.tab", delimiter="\t")

df.columns = ['Chr', 'Start', 'ExonicFunc', 'Exon', 'Protein_Change']

list_of_pros = df['Protein_Change'].values.tolist()

conversion = convert_gba_notation(list_of_pros)

conversion.to_csv(r'notation_conversion.tab', index = False, sep="\t")


