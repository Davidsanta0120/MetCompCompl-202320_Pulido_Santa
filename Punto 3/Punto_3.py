import numpy as np
import sympy as sp

x=sp.Symbol("x")

def GetLegendreRecursive(n,x):

    if n==0:
        poly = sp.Number(1)
    elif n==1:
        poly = 1-x
    else:
        poly =((((2*(n-1))+1-x)*GetLegendreRecursive(n-1,x))-((n-1)*GetLegendreRecursive(n-2,x)))/n
    return sp.expand(poly,x)
print(GetLegendreRecursive(5,x))
def GetNewton(f,df,xn,itmax=10000,precision=1e-14):
    
    error = 1.
    it = 0
    
    while error >= precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn)
            
            error = np.abs(f(xn)/df(xn))
            
        except ZeroDivisionError:
            print('Zero Division')
            
        xn = xn1
        it += 1
        
    if it == itmax:
        return False
    else:
        return xn
    
def GetDLegendre(n,x):
    Pn = GetLegendreRecursive(n,x)
    return sp.diff(Pn,x,1)
    
def GetRoots(f,df,x,tolerancia = 10):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)

        if  type(root)!=bool:
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
    return Roots
def GetAllRootsGLeg(n):
    j=n+((n-1)*np.sqrt(n))
    xn = np.linspace(-1,j,100)
    
    Legendre = []
    DLegendre = []
    
    for i in range(n+1):
        Legendre.append(GetLegendreRecursive(i,x))
        DLegendre.append(GetDLegendre(i,x))
    
    poly = sp.lambdify([x],Legendre[n],'numpy')
    Dpoly = sp.lambdify([x],DLegendre[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots

print(GetAllRootsGLeg(6))
def GetWeightsGLag(n):
    l_n_1 = GetLegendreRecursive(n+1,x)
    l_n_1 = sp.lambdify(x,l_n_1)
    x_k = GetAllRootsGLeg(n)
    c_k = np.array([])
    for i in x_k:
        p = (i)/(((n+1)**2)*((l_n_1(i))**2))
        c_k = np.append(p,c_k)
    return c_k
print(GetWeightsGLag(3))

def GetHermitRecursive(n,x):

    if n==0:
        poly = sp.Number(1)
    elif n==1:
        poly = 2*x
    else:
        poly = (2*x*(GetHermitRecursive(n-1,x)))-(2*(n-1)*(GetHermitRecursive(n-2,x)))
    return sp.expand(poly,x)
print(GetHermitRecursive(7,x))
def GetDHermit(n,x):
    Pn = GetHermitRecursive(n,x)
    return sp.diff(Pn,x,1)
def GetAllRootsHermit(n):
    j=np.sqrt(4*n+1)
    xn = np.linspace(-j,j,100)
    
    Hermit = []
    DHermit = []
    
    for i in range(n+1):
        Hermit.append(GetHermitRecursive(i,x))
        DHermit.append(GetDHermit(i,x))
    
    poly = sp.lambdify([x],Hermit[n],'numpy')
    Dpoly = sp.lambdify([x],DHermit[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots
print(GetAllRootsHermit(4))