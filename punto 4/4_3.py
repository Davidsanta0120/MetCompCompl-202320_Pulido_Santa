import numpy as np 
def problema4_3 ()->float:
    a=-1
    b=1
    s=100000
    x = np.random.uniform(a,b,s)
    y = np.random.uniform(a,b,s)
    z = np.random.uniform(a,b,s)
    
    omega=np.sqrt(x**2+y**2+z**2)
    
    
    f= lambda Om: np.sin(Om**2)*np.exp(Om**2)
    summ=0
    for i in omega:
        if i <=1:
            fs=f(i)
            summ+=fs
    
    
    I=(b-a)/s*summ
    return I

print(problema4_3())
