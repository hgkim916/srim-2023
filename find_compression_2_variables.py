from fpylll import *
from fpylll.fplll.gso import MatGSO
import math

def LatticeBasis2Variables(degree,box_x,box_y):
    number_of_rows = math.comb(degree + 2, 2)
    number_of_columns = box_x*box_y
    basis = IntegerMatrix(number_of_rows,number_of_columns)
    i = 0
    for current_degree in range(degree+1):
        for x_bot in range(current_degree+1):
            y_bot = current_degree - x_bot
            j = 0
            for x_top in range(box_x):
                for y_top in range(box_y):
                    #print(x_top,x_bot,y_top,y_bot)
                    basis[i,j] = math.comb(x_top,x_bot)*math.comb(y_top,y_bot)
                    j += 1
            i += 1
    return basis

def PossibleDynCompressing2Variables(degree,box_size,number_of_solutions,
                                     show_only_dynamically_compressing = True,
                                     show_vector = False,
                                     force_double_dependence = False):
    # Input: degree, box_size, #of solutions
    # Output: coordinates of each polynomial which can possibly send [box_size]^2 to [box_size]. 
    # note not all of these vectors will work!
    number_of_rows = math.comb(degree + 2, 2)
    number_of_columns = box_size*box_size
    if box_size<=degree:
        print("invalid box size!")
        return

    A = LatticeBasis2Variables(degree,box_size,box_size)
    
    M = MatGSO(A)
    _ = M.update_gso()
    enum = Enumeration(M,strategy = EvaluatorStrategy.BEST_N_SOLUTIONS,nr_solutions = number_of_solutions,sub_solutions=False)
    euclidean_upper_bound_squared = (number_of_columns)*(math.floor((box_size)/2)**2)
    #print(euclidean_upper_bound_squared)
    e1 = enum.enumerate(0, number_of_rows, euclidean_upper_bound_squared, 0)
    count = 0
    found_vectors = []
    print(len(e1))
    for i,j in e1:
        count += 1
        if count%10000 == 0:
            print(count)
        # Check for unwanted cases
        if max(j[3:]) == 0 and min(j[3:]) == 0:
            print(count,"(linear solution)")
            continue
        if max(j[-(degree+1):]) == 0 and min(j[-(degree+1):]) == 0:
            #print(count,"(lower degree solution)")
            continue
        #print((i),tuple(map(int,j)))
        
        # Check if we found the same polynomial up to a constant
        non_constant_part = j[1:]
        if non_constant_part in found_vectors:
            #print(count,"(already found up to a constant)")
            continue
        else:
            found_vectors.append(non_constant_part)
        
        # Find actual vector
        x_dependence = False
        y_dependence = False
        vector = IntegerMatrix(1,number_of_columns)
        coord = 0
        for current_degree in range(degree+1):
            for x_bot in range(current_degree+1):
                y_bot = current_degree - x_bot
                if int(j[coord]) != 0:
                    if x_bot != 0:
                        x_dependence = True
                    if y_bot != 0:
                        y_dependence = True
                vector[0].addmul(A[coord],int(j[coord]))
                coord += 1
            #print(A[coord],int(j[coord]))
            #print(vector)
            
        # Turn into a Python list
        new_vector = [0]*(number_of_columns)
        for coord in range(number_of_columns):
            new_vector[coord] = vector[0,coord]
        #print(new_vector)
        
        if show_only_dynamically_compressing and max(new_vector)-min(new_vector)+1 > box_size:
            #print(count,"(doesn't dynamically compress)")
            continue
        
        if force_double_dependence and (not x_dependence or not y_dependence):
            #print(count,"(doesn't depend on both x and y)")
            continue
            
        print(count,i,"\tCompression:",box_size,"to",max(new_vector)-min(new_vector)+1,"\tCoordinates:",tuple(map(int,j)))#,"\tVector:",new_vector)

print("Hello World!") # So you know the program is running
#A = LatticeBasis2Variables(2,8,8)
#print(A)
#LLL.reduction(A)
#print(A)
PossibleDynCompressing2Variables(13,19,1000,force_double_dependence = True,show_vector=False)