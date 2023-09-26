import numpy as np 
import matplotlib.pyplot as plt

def problema4_1 ():
    a=0
    b=np.pi
    
    N = list(np.linspace(1000,100000,150))
    
    f=lambda x: np.exp(-x)*np.sin(x)
    valreal=0.5+np.exp(-np.pi)/2
    error_N=[]
    
    for i in range(0,len(N)): 
        s=int(N[i])
        x = np.random.uniform(a,b,s)
        fs=f(x)
        I=(b-a)/s*np.sum(fs)
        error= abs(valreal-I)/valreal*100
        error_N.append(error)
    
    plt.scatter(N,error_N)
    plt.show()

problema4_1()