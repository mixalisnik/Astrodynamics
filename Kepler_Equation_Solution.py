#Initially, import all neccessary libraries.
import math as m
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def lagrange_t1(e,nmax): #This function takes input eccentricity e and max order nmax,
                         #that will be used to compute the series expansion.
  x=sp.Symbol('x')       #It implements Lagrange's Expansion Theorem, giving output of y(x) (i.e. E(M)) accordingly.
  y=x
  alpha=e
  phi=sp.sin(x)
  
  for n in range(1,nmax+1): #Series Calculation upto nmax order
    phin=phi**n
    temp=alpha**n/m.factorial(n)*(phin.diff(x,n-1))
    y=y+temp
    
  
  return y #y is a symbolic expression of x

def plotlgr(e,nmax): #This function plots the solution of Kepler's Equation using the Lagrange Theorem
                     #for eccentricity e and max order nmax.
    
    f1=lagrange_t1(e,nmax)  #f1 (i.e y) is a symbolic expression. To use it as a callable function, we use lambdify. 
    f11=sp.lambdify('x',f1,'numpy') 
    x1=np.linspace(0,2*np.pi,1000) #Create an equally spaced vector of x in range [0,2π], containing 1000 elements.
    a_matrix=np.zeros([np.size(x1),2]) #Initialize a_matrix of 1000x2 dimensions. 1st column will store x (i.e. M values), 
                                       #2nd column will store y (i.e. E(M) values)
    j1=0                                
    for i in x1:  #For each element in x1 vector, find its corresponding y (f11(x)) and store x and y in a_matrix
        a_matrix[j1,0]=i #x Values
        a_matrix[j1,1]=f11(i) #y Values
        j1=j1+1
        
    plt.plot(a_matrix[:,0],a_matrix[:,1],label=('e='+str(e)+' n='+str(nmax))) #Plot (M,E(M)) with corresponding labels
    plt.legend() #Display legend
    plt.title('Solution of Kepler\'s Equation using Lagrange\'s Theorem ') #Set title of plot.
    plt.ylabel(' E(M) ') #Set x and y labels
    plt.xlabel('M')
    return 

def jbessel(n,kmax): #This function computes the Bessel function of order n, keeping terms up to kmax. 
                      #The output of this function (y) (i.e. Series Expansion Result up to kmax terms) is a symbolic expression of x.
  x=sp.Symbol('x')
  y=0
  
  for k in range(0,kmax):
    
    temp=(-1)**k/(m.factorial(k)*m.factorial(n+k))*(x/2)**(n+2*k)
    y=y+temp
    
  
  return y

def plotbessel(n,kmax): #Plots the Bessel functions of order n, keeping terms up to kmax.
    
    f1=jbessel(n,kmax) #f1 is a symbolic expression of x, that needs to be converted to a callable function of x using lambdify.
    f11=sp.lambdify('x',f1,'numpy')
    x1=np.linspace(0,15,1000) #Create an equally spaced vector of x in range [0,2π], containing 1000 elements.
    a_matrix=np.zeros([np.size(x1),2])  #Initialize a_matrix of 1000x2 dimensions. 1st column will store x, 
                                       #2nd column will store y (Jn(x))
    j1=0
    for i in x1:  #For each element in x1 vector, find its corresponding y (f11(x)) and store x and y in a_matrix
        a_matrix[j1,0]=i #x
        a_matrix[j1,1]=f11(i) #y
        j1=j1+1
        
    plt.plot(a_matrix[:,0],a_matrix[:,1],label='n='+str(n)) #Plot (x,Jn(x)) with corresponding labels
    plt.title('Bessel Functions Jn')
    plt.ylabel(' Jn(x) ') #Axes labels
    plt.xlabel('x')
    plt.legend() #Display legend
    plt.grid() #Display grid
    return 

def fourier_bessel_s1(kmax,jmax,e): #This function computes the symbolic expression of the Fourier - Bessel Series Expansion
                                   #for eccentricity e, specifying j-max and k-max.
    x=sp.Symbol('x')
    y=x
    for j in range(1,jmax): 
      f1=jbessel(j,kmax) #Find the Bessel function corresponding to j order, keeping upto kmax terms.
      f11=sp.lambdify('x',f1,'numpy') #Use lambdify to convert f1 to a callable function of x.
      temp=2/j*f11(j*e)*sp.sin(j*x) #Implement the Fourier Bessel Series Expansion
      y=y+temp
    return y #y is a symbolic expression of x

def plot_fb(kmax,jmax,e): #Plots the Fourier-Bessel function
    
    f1=fourier_bessel_s1(kmax, jmax, e)
    f11=sp.lambdify('x',f1,'numpy')
    x1=np.linspace(0,2*np.pi,1000)
    a_matrix=np.zeros([np.size(x1),2])
    j1=0
    for i in x1: 
        a_matrix[j1,0]=i #x
        a_matrix[j1,1]=f11(i) #y
        j1=j1+1
        
    plt.plot(a_matrix[:,0],a_matrix[:,1],label=('j='+str(jmax)+' e='+str(e)+' k='+str(kmax))) #Plot with corresponding labels
    plt.legend() #Display legend
    plt.title('E(M) using Fourier - Bessel Expansion') #Set title
    plt.ylabel(' E(M) ') #Set axes' labels
    plt.xlabel('M')
    return 



def f(x,mvi,e): #Function that will be used as f(x) for Newton Rhapson Method.
    return x-e*np.sin(x)-mvi

def dfdx(x,e): #Derivative of f(x) with respect to x
    return 1-e*np.cos(x)
def kepler_newton(e1,acc,xi1): #This function computes and plots the Solution of Keplers Equation, using Newton Rhapson Method.
    #Takes input matrix of eccentricities e1, desired accuracy of computed root (acc), initial root guess (xi1)
    mv=np.linspace(0,2*np.pi,1000) #We need to find the roots of f(x) for each M in [0,2π]
    
    n1=np.size(mv)
    
    nmax=15 #Maximum number of iterations that guards from slow convergence.
    a=np.zeros((n1,2)) 
    
    for k in e1: #e1 is regarded as a vector of eccentricities. Thus, kepler_newton function can receive multiple values of e at once.
      j1=0 #Attention: Since e1 is regarded as a matrix, even a single value of e must be in the form of [e]. 
      xi=xi1 #Initial guess of root of f(x). Neccessary for Newton Rhapson initialization
      for j in mv: 
          
      
          for i in range(1,nmax+1): #Implementation of Newton Rhapson Method
              
              if abs(dfdx(xi,k))<10**-6: #If derivative is almost zero, error is displayed and for loop is terminated, displaying following message
                  print('Error: Derivative at xi='+str(xi)+' close to zero. df/dx at xi equals '+str(dfdx(xi)))
                  break
              if i==nmax+1: #If max iterations are reached without reaching desired accuracy, error message is displayed.
                  print('Error: Maximum number of iterations reached without reaching desired accuracy.')
                  
              xi1=xi-f(xi,j,k)/dfdx(xi,k) #Implementation of Newton Rhapson Method
              
              if abs(xi1-xi)<acc: #Termination criterion of Method. 
                  a[j1,0]=xi1 # E value (x) for M value (mv)
                  a[j1,1]=j # M value
                  break
              
              xi=xi1 #Update of xi value
              
          j1=j1+1
      y2=a[:,0];
      x2=a[:,1];
      plt.plot(x2,y2,label='e='+str(k)) #Plot (M,E(M))
      
    
    temp=np.linspace(0,2*np.pi,13) #Plot Configurations
    plt.legend()
    plt.title('Solution of Kepler\'s Equation using Newton - Rhapson Method')
    plt.ylabel('E values ')
    plt.xlabel('M values')
    plt.xticks(temp,['0','π/6','π/3','π/2','2π/3','5π/6','π','7π/6','4π/3','3π/2','5π/3','11π/6','2π'])
    plt.yticks(temp,['0','π/6','π/3','π/2','2π/3','5π/6','π','7π/6','4π/3','3π/2','5π/3','11π/6','2π'])
    plt.grid()
    
    return 

#Task 1  - Newton Rhapson
plt.figure(0) 
kepler_newton([0.1,0.3,0.5,0.7,0.9],10**-4,2)

#Task 2 - Lagrange Theorem
plt.figure(1) #e=0.3
kepler_newton([0.3],10**-4,2)
plotlgr(0.3,3)
plotlgr(0.3,10)
#for e=0.3, Lagrange's Theorem gives a very precise solution of Kepler's Equation. Both for n=3 and n=10,
# this approach converges extremely well with Newton - Rhapson Method for e<<.
plt.figure(2) #e=0.9
kepler_newton([0.9],10**-4,2)
plotlgr(0.9,3)
plotlgr(0.9,10)
#for e=0.9, Lagrange's Theorem begins to diverge from Newton - Rhapson results. 
#Both for n=3 and n=10, the Error of this method is smaller around Pi, increasing as we move further away from Pi.
#Lagrange's Theorem demonstrates a periodicity of the Solution, and is therefore not accurate and precise enough for e>>.
#As n increases, the period of the solution becomes smaller, resulting in more frequent maximum and minimum values.

#Task 3 - Fourier Bessel Expansion

plt.figure(3) #Bessel Functions Plot Notice:Functions are called for kmax=100. Lower value as desired to avoid overflows.
for i1 in range(1,6):
    plotbessel(i1,100)

plt.figure(4) #e=0.3 Notice:Functions are called for kmax=100. Lower value as desired to avoid overflows.
kepler_newton([0.3],10**-4,2)
plot_fb(100,3,0.3)
plot_fb(100,10,0.3)

plt.figure(5) #e=0.9 Notice:Functions are called for kmax=100. Lower value as desired to avoid overflows.
kepler_newton([0.9],10**-4,2)
plot_fb(100,3,0.9)
plot_fb(100,10,0.9)
#for e=0.3, Fourier - Bessel Expansion gives a very precise solution of Kepler's Equation. 
#Both for n=3 and n=10, this approach converges extremely well with Newton - Rhapson Method,
#as well as with Lagrange's Theorem.
#for e=0.9, Fourier - Bessel Expansion begins to diverge from Newton - Rhapson results. 
#For n=3 the error is significant and demonstrates periodicitiy, although the approach is more accurate than Lagrange's Theorem.
#For n=10 the error is infinitesimal, and the approach converges with the solution of Newton - Rhapson Method.
#To conclude, Newton Rhapson method and Fourier Bessel Expansion for n=10 seem to yield the most precise results, and therefore must be preferred over Lagrange's Theorem.
plt.show()


