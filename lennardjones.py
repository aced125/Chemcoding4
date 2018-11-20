#Optimisation of n particles in Lennard-Jones potential with Covariance Matrix Adapted-Evolution Strategy (CMA-ES)



import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from operator import add
import time
import cma
#import cma.test


#Testing cma module
#print('\nTesting CMA module; to read more, check out https://arxiv.org/abs/1604.00772\n')
#cma.test.main()

#Initial conditions

x0 = np.array((0,0,0))
x1 = np.array((1,0,0))
x2 = np.array((0,1,0))
x3 = np.array((0,0,0.9))
x4 = np.array((1,1,0))
x5 = np.array((1,0,1))
x6 = np.array((0,1,1))


particle_number = int(input('\nHow many particles is your system?\n'))

iterations_to_search_for = int(input('\nHow many searches would you like to conduct? (Recommended 15)\n'))



gradstep = -10**-3

coordlist = [x0,x1,x2,x3,x4,x5,x6]

def ld(r1,r2):
    return linalg.norm(r1-r2)

def LJpotential(r1,r2):
    return 4*(ld(r1,r2)**-12 - ld(r1,r2)**-6)


def sumofLJpotential(coordlist):
    
    sum = 0
    for i in range(len(coordlist)):
        for j in range(i+1,len(coordlist)):
            sum+=LJpotential(coordlist[i],coordlist[j])
            
    return sum


def LJgrad_wrt_rk(k):
    
    sum = np.array((0,0,0),dtype = float)
    for j in range(len(coordlist)):
        if j != k:
            sum += (coordlist[k]-coordlist[j])*(-2*ld(coordlist[k],coordlist[j])**-14 + ld(coordlist[k],coordlist[j])**-7)
    
    
    return 24*sum



#gradlist = [LJgrad_wrt_rk(k) for k in range(len(coordlist))]


#iteration_number = 0

#print('\nOptimizing geometry, please wait...\n')

#start = time.time()

#for iteration in range(10**4):
#while linalg.norm(gradlist) > 10**-10:
    
    #gradlist = [LJgrad_wrt_rk(k) for k in range(len(coordlist))]

    #smallgrad = list(map(lambda x:x*gradstep,gradlist))

    #coordlist = list(map(add,coordlist, smallgrad))
    
    #if iteration_number%1000==0:
    
        #print(f'iterations:{iteration_number}')
        #print(f'norm of gradient vector: {linalg.norm(gradlist)}')
    #iteration_number+=1

#end = time.time()

#print(f'\nIterations required for convergence: {iteration_number}')
#print(f'\nTime taken: {round(end-start,4)} seconds')
#print(f'\nTime for one iteration: {(end-start)/iteration_number} seconds')
#print(f'\nconverged potential: {round(sumofLJpotential(coordlist),5)}')
    

def list_to_coord_conversion(listcoordlist):
    """Converts list of n coordinates in R3 to n/3 vector coordinates"""

    arraycoordlist = [np.array(listcoordlist[3*i:3*(i+1)]) for i in range(particle_number)] 
    return arraycoordlist

def coord_to_list_conversion(coordlist):
    """Converts n/3 vector coordinates in R3 to n list coordinates"""

    listcoordlist = [coordlist[i][j] for i in range(len(coordlist)) for j in range(len(coordlist[0]))]
    return listcoordlist


def LJpotentialsum(list1):
    """Computes the Lennard-Jones potential given a list of n coordinates"""
    arraycoordlist = [np.array(list1[3*i:3*(i+1)]) for i in range(int(len(list1)/3))]
    return sumofLJpotential(arraycoordlist)


xoptlist = []
minpotentiallist = []
for i in range(iterations_to_search_for):
    
    xopt, es = cma.fmin2(LJpotentialsum, list(np.random.randn(particle_number*3)), 0.5)
    xoptlist.append(xopt)
    minpotentiallist.append(es.result[1])
    print(f'\niteration: {i}, which yielded a potential of: {es.result[1]} \n ')


min_index = np.argmin(np.array(minpotentiallist))
finalcoordlist = list_to_coord_conversion(xoptlist[min_index])


def origintranslator(finalcoordlist):
    """"Centres the solution so that one of the points is the origin"""
    
    centre = finalcoordlist[0]
    centredfinalcoordlist = finalcoordlist - centre
    return centredfinalcoordlist


finalcoordlist = list_to_coord_conversion(xopt)
centredfinalcoordlist = origintranslator(finalcoordlist)
coordmatrix = np.array(centredfinalcoordlist)


print(f'\nThe minimum potential CMA-ES was able to find was: {min(minpotentiallist)}\n in units of epsilon')

print('\nWriting an xyz file to the current directory\n')


fn = 'LJpotential'+str(particle_number)+'particles.xyz'
file = open(fn,'w')
file.write(f'{particle_number}\n')
file.write(f'Optimum geometry of {particle_number} particles minimised by Lennard Jones potential\n')

for i in range(particle_number):
    file.write(f'particle{i} {coordmatrix[i][0]} {coordmatrix[i][1]}  {coordmatrix[i][2]}\n')

file.close()


print('\nA plot of these points is shown here\n')



#Plot position of particles in equilibrium geometry

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs = coordmatrix[:,0]
ys = coordmatrix[:,1]
zs = coordmatrix[:,2]
ax.scatter(xs, ys, zs, c='r', marker='o')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.title('Position of particles')

plt.show()