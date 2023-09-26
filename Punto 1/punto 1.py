import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
f = lambda t: (np.pi)*((0.25/2)**2)*(0.05)*(np.cos(3.5*t))*(np.cos(2*(np.pi)*(1000)*t))
#z = lambda t: (-1/1750)*((np.pi)*((0.25/2)**2)*0.05*(((-3.5)*(np.sin(3.5*t))*(np.cos(2*t*1000*np.pi)))-((2*1000*np.pi)*(np.sin(2*t*1000*np.pi))*(np.cos(3.5*t)))))
z = lambda t: (-1/1750) * (Deriv_flujo_central(t,f))
x = np.linspace(-0.0001,0.1,num=5000)
def Deriv_flujo_central(x,f,h=1e-6):
    return (f(x+h)-f(x-h))/(2*h)

def GetNewtonMethod(f,df,xn,itmax=100,precision=1e-8):
    
    error = 1.
    it = 0
    
    while error > precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn,f)
            # Criterio de parada
            error = np.abs(f(xn)/df(xn,f))
            
        except ZeroDivisionError:
            print('Division por cero')
            
        xn = xn1
        it += 1
        
   # print('Raiz',xn,it)
    
    if it == itmax:
        return False
    else:
        return xn
def GetAllRoots(x, tolerancia=10):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewtonMethod(z,Deriv_flujo_central,i)
        
        if root != False:
            
            croot = np.round(root, tolerancia)
            
            if (croot not in Roots)and(croot>=0):
                Roots = np.append(Roots,croot)
                
                
    Roots.sort()
    
    return Roots

I=z(x)
x_1 = np.arange(0.0001,5)
roots = (GetAllRoots(x))[0:4]
plt.plot(x[:100],I[:100])
plt.scatter(roots, z(roots), c='red')
plt.show()
print("Los valores en donde la corriente se hace 0 son: ", roots )
