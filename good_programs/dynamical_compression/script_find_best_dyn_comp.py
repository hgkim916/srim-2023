from find_compression import *

# Tries to find polynomials of degree d which compress [n] for increasing n until 4 consecutive failures.

### Change variables as desired

# Range of d
d_min = 2
d_max = 14

consecutive_failures_limit = 4 # How many consecutive failures needed before program stops checking

###

results = []
for d in range(d_min,d_max+1):
    last_success = 6
    compression = 7
    while compression-last_success <= consecutive_failures_limit:
        succeeded = print_dynamical_compression(d,compression,1000000,show_only_dyn_comp = True,show_vector = False,print_failures=False)
        if succeeded:
            last_success = compression
        compression += 1
    print("best n:",last_success+d)
    results.append(last_success+d)

print("----------")
for d in range(d_min,d_max+1):
    print("d:",d,"\tbest n:",results[d-2])