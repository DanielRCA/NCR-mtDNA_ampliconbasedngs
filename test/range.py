

import sys
import csv
from itertools import groupby
from operator import itemgetter


data=sys.argv[1]
the_input = []
with open (data) as rfile:
    for position in rfile:
        the_input.append(int(position.strip()))
rfile.close()
    

def split_list(n):

    return [(x+1) for x,y in zip(n, n[1:]) if y-x != 1]

def get_sub_list(my_list):

    my_index = split_list(my_list)
    output = list()
    prev = 0
    for index in my_index:
        new_list = [ x for x in my_list[prev:] if x < index]
        output.append(new_list)
        prev += len(new_list)
    output.append([ x for x in my_list[prev:]])
    return output

full_list = []
for group in get_sub_list(the_input):
    if group[0] == group [-1]-1:
        x = (group[0],";",group[-1],";")
        full_list.append(x)
    elif group[0] == group [-1]:
        x = (group[0],";")
        full_list.append(x)
    else:
        x = (group[0],"-",group[-1],";")
        full_list.append(x)
        
print(full_list)



        
