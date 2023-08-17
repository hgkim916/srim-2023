from find_compression import *

for d in range(2,15):
    last_success = 6
    compression = 7
    while compression-last_success <= 4:
        succeeded = print_dynamical_compression(d,compression,1000000,show_only_dyn_comp = True,show_vector = True)
        if succeeded:
            last_success = compression
        compression += 1
    print("best n:",last_success+d)
