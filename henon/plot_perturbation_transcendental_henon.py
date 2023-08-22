import numpy as np
import matplotlib.pyplot as plt

# list of coefficients for polynomials s_d for d=0,1,...,9
# (where zeros are skipped and leading term has positive coeff)
SD_LIST = [
    [1],
    [1],
    [-9/8,1/2],
    [-7/6,1/6],
    [147/128,-29/48,1/24],
    [6/5,-5/24,1/120],
    [-1181/1024,7219/11520,-31/576,1/720],
    [-169/140,157/720,-1/90,1/5040],
    [37827/32768,-45277/71680,2621/46080,-11/5760,1/40320],
    [1523/1260,-19987/90720,41/3456,-17/60480,1/362880]
] 

###

def s(x):
    '''
    computes s(y), where s is pointwise limit of s_d, namely s(y)=2/sqrt(3)*sin(pi*y/3)
    '''
    if type(x) is int:
        if (x %6)%3==0:
            return 0
        elif (x % 6)-3 < 0:
            return 1
        else: return -1
    else:
        return 2/np.sqrt(3) * np.sin(np.pi*x/3)


def sd(y,d):
    '''
    computes polynomial s_d(y) for d <= 9
    '''
    if d == np.Inf:
          return s(y)
    elif d > 9:
        return ValueError # not implemented yet
    coeffs = SD_LIST[d]
    flip = -1 if d%4>1 else 1
    return sum([c*pow(y,2*k) for k,c in enumerate(coeffs)])*(y**(d%2))*flip


def f(x,y,d=np.Inf):
    '''
    Henon map f_d(x,y)=(y,-x+s_d(y))
    '''
    return (y,-x + sd(y,d))

###

def draw_ellipse(centre,pass_through,mu): 
	'''
	draws an ellipse centred at (centre), passing through (pass_through)
	which arises from forward orbit of perturbation near periodic point with Jacobian [[0,1],[-1,(mu)]]
	'''
	x0,y0 = centre
	x,y = pass_through
	C = ( (x-x0) + mu*(y-y0))**2/(1-mu**2) + (y-y0)**2
	xlist, ylist = [], []
	for i in range(-999,1000):
		yi = y0 + np.sqrt(C) * i/1000
		xi = x0 -mu*(yi-y0) + ((1-mu**2)*(C-(yi-y0)**2))**0.5
		xlist.append(xi)
		ylist.append(yi)
		xi = x0 -mu*(yi-y0) - ((1-mu**2)*(C-(yi-y0)**2))**0.5
		xlist.append(xi)
		ylist.append(yi)
		
	plt.scatter(xlist,ylist)
	plt.plot()
	plt.show()
        

###

def find_period(x,y,d=np.Inf):
    '''
    finds period of integer point (x,y) for f_d
    '''
    cnt = 0
    z,w = x,y
    while 1:
        cnt += 1
        z,w = f(z,w,d)
        print(z,w)
        if (z,w) == (x,y):
            return cnt
        if cnt >= 1e+5:
            raise ValueError # not a periodic point


###



def iteration_plot(x,y,n=1,d=np.Inf,eps=10000,plot_from=0,color=False,dot_size=40):
    '''
    plots forward orbit of point (x,y) under the Henon map f_d
    :param n: iterates function (f_d)^n
    :param eps: number of forward iterations
    :param plot_from: plot points from (plot_from)th iterate to (eps)th iterate
    :param color: color of dot, specified by colors in matplotlib
    '''
    x_list, y_list = [],[]
    for i in range(eps):
        if i >= plot_from:
            x_list.append(x)
            y_list.append(y)
        for _ in range(n):
            x,y = f(x,y,d)
    if color is False:
        plt.scatter(x_list,y_list,s=dot_size)
    else:
        plt.scatter(x_list,y_list,s=dot_size,c=color)
    plt.plot()



def plot_perturbations_from_transcendental():
    '''
    plots forward orbits of small perturbations from integer points near the origin, under iteration of the transcendental Henon map
    '''
    iteration_plot(0,0.15,eps=200000,color='r',dot_size=2)
    iteration_plot(1,0.1,eps=200000,color='b',dot_size=2)
    iteration_plot(-0.99,1,eps=200000,color='m',dot_size=2)
    iteration_plot(0,2.02,eps=200000,color='g',dot_size=2)
    iteration_plot(3,0.1,eps=200000,color='k',dot_size=2)
    iteration_plot(2,3.07,eps=200000,color='c',dot_size=2)
    iteration_plot(4,3.07,eps=200000,color='darkblue',dot_size=2)
    iteration_plot(3,3.15,eps=200000,color='crimson',dot_size=2)
    iteration_plot(1.1,4.05,eps=200000,color='darkviolet',dot_size=2)
    iteration_plot(1,4.05,eps=200000,color='lawngreen',dot_size=2)

    plt.grid(visible=True)
    plt.rcParams.update({'font.size': 15})
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.title("Perturbations near integer points")

    N = 7
    plt.xticks(np.arange(-N-3,N+4))
    plt.yticks(np.arange(-N-3,N+4))
    plt.xlim(-N,N)
    plt.ylim(-N,N)

    plt.show()