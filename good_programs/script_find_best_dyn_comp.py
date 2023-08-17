from find_compression import *

# Tries to find polynomials of degree d which compress [n] for increasing n until 4 consecutive failures.

results = []
for d in range(2,15):
    last_success = 6
    compression = 7
    while compression-last_success <= 4:
        succeeded = print_dynamical_compression(d,compression,1000000,show_only_dyn_comp = True,show_vector = False,print_failures=False)
        if succeeded:
            last_success = compression
        compression += 1
    print("best n:",last_success+d)
    results.append(last_success+d)

print("----------")
for d in range(2,15):
    print("d:",d,"\tbest n:",results[d-2])