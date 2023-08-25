# Demonstrates the next_list_colex_order function.

from henon_search_for_longer_cycles import next_list_colex_order

list_length = 4
min_entry = -2
max_entry = 2
a = [min_entry]*4

while a != None:
    print(a)
    a = next_list_colex_order(a,min_entry,max_entry)
