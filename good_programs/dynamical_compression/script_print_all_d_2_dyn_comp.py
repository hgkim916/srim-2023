from find_compression import *

for compression in range(1,7):
    print_dynamical_compression(2,compression,2000,show_vector=False,show_translated_poly=True,print_failures=False)