from math import gcd, ceil
import itertools
from scipy import sparse
import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon

#Calculate equilibrium matrix B
def calcB(Nd, Cn, dof):
    m, n1, n2 = len(Cn), Cn[:,0].astype(int), Cn[:,1].astype(int)   #m is number of elements, n1 and n2 are list of nodes
    l, X, Y = Cn[:,2], Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1]         #l is the length (given in Cn), X and Y are X and Y run
    d0, d1, d2, d3 = dof[n1*2], dof[n1*2+1], dof[n2*2], dof[n2*2+1] #Finding the 4 forces (n1x, n1y, n2x, n2y) correspo to the each element
    s = np.concatenate((-X/l * d0, -Y/l * d1, X/l * d2, Y/l * d3))  #force on the nodes caused by 1 unit of tension in the bar
    r = np.concatenate((n1*2, n1*2+1, n2*2, n2*2+1))                #List of indexes in the B matrix to store the forces
    c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m))) #Creates a list from 1 to number of elements
    return sparse.coo_matrix((s, (r, c)), shape = (len(Nd)*2, m))  #assembling s values into the B matrix

#Solve linear programming problem
def solveLP(Nd, Cn, f, dof, st, sc, jc):
    l = [col[2] + jc for col in Cn]             #Setting length, adjusting by jc. Joint cost is the cost of joints!
    B = calcB(Nd, Cn, dof)                      #Solving for the B matrix
    a = cvx.Variable(len(Cn))                   #Variable which expresses member cross-sectional areas
    obj = cvx.Minimize(np.transpose(l) @ a)     #Objective function
    q, eqn, cons= [],  [], [a>=0]               #q is the force vector, eqn will be force-equilibrium equation, const will be where the constraints go for the optimization problem
    if case == 2:                               #Adding max area constrain for case 2
       cons.extend([a<=.5])
    for k, fk in enumerate(f):                  #Looping through the diferent load cases
        q.append(cvx.Variable(len(Cn)))         #Adding the force for each member as an optimization variable
        eqn.append(B @ q[k] == cvx.multiply(fk,dof))        #Main force balance equation
        cons.extend([eqn[k], q[k] >= -sc * a, q[k] <= st * a])  #Adding the optimization constraints
    prob = cvx.Problem(obj, cons)               #Defining the optimization problem. Default solver to be ECOS
    vol = prob.solve()                          #Solving porblem
    q = [np.array(qi.value).flatten() for qi in q] #write member forces
    a = np.array(a.value).flatten()             #Write section areas
    return vol, a, q

#Visualize truss
def plotTruss(Nd, Cn, a, q, threshold, string, update = True):  #Literally taken directly from the paper
    plt.ion() if update else plt.ioff()
    plt.clf(); plt.axis('off'); plt.axis('equal'), plt.draw()
    plt.title(string)
    tk = 5 / max(a)
    for i in [i for i in range(len(a)) if a[i] >= threshold]:
        if all([q[lc][i]>=0 for lc in range(len(q))]): c = 'b'
        elif all ([q[lc][i] <=0 for lc in range(len(q))]) : c = 'r'
        else: c='tab:gray'
        pos = Nd[Cn[i, [0,1]].astype(int),:]
        plt.plot(pos[:,0], pos[:,1], c, linewidth = a[i]*tk)
        if text == 1:
            [xavg, yavg] = (pos[1,:]+pos[0,:])*.5; plt.text(xavg, yavg, "F="+str(round(q[0][i],3))+"\nA="+str(round(a[i],3)),weight='extra bold')
    plt.pause(.01) if update else plt.show()

#Visualize possible grid layout
def plotTrussGrid(Nd, Cn):
    plt.ion()
    plt.clf(); plt.axis('off'); plt.axis('equal'), plt.draw()
    plt.title("Possible Member Layout, Nodes, DOFs")
    tk = 5
    for i in range(int(np.size(Cn) / 4)): #Looping over the members
        c='tab:gray'
        pos = Nd[Cn[i, [0,1]].astype(int),:] #Get coordinates of both nodes
        plt.plot(pos[:,0], pos[:,1], c, linewidth = 1)
        if text == 1:
            [xavg, yavg] = (pos[1,:]+pos[0,:])*.5; plt.text(xavg, yavg, "b"+str(i+1),color='g')
    if text == 1:
        for i in range(int(np.size(Nd)/2)):
            nodecor = Nd[i, [0,1]]
            plt.text(nodecor[0], nodecor[1], "n"+str(i+1),color='r',weight='extra bold')
            plt.text(nodecor[0], nodecor[1], "dof" + str((i + 1)*2-1)+"   ", color='b',horizontalalignment='right',verticalalignment='center')
            plt.text(nodecor[0], nodecor[1], "\ndof" + str((i + 1)*2), color='b',verticalalignment='top', horizontalalignment='center')
    plt.show()
    plt.pause(4)

def trussopt(width, height, st, sc, jc, string):
    if case == 1 or case == 2 or case == 3:
        poly = Polygon([(0,0), (width, 0), (width, height), (0, height)]) #Create the area over which to iterate
    elif case == 4:
        poly = Polygon([(0, 0), (width/4, 0), (width/4, height/2), (3*width/4, height/2), (3*width/4, 0), (width,0), (width, height), (0, height)])
    convex = True if poly.convex_hull.area == poly.area else False  #checks weather area is convex
    xv, yv = np.meshgrid(range(width+1), range(height+1))  #defines a meshgrid with the possible points
    pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)] #crete a list of points for the possible nodes
    Nd = np.array([[pt.x,pt.y] for pt in pts if poly.intersects(pt)]) #Creating a list of the possible nodes (note the the "if" statement excludes points that would be within a central nogo zone
    dof, f, PML = np.ones((len(Nd), 2)), [], [] #Setting a dof vector to be an array of ones, NX2, also defining f vec and PML vec (not sure what these are?)
    #Load and support conditions
    for i, nd in enumerate(Nd): #Looping through all the nodes one by one
        if case == 4:
            if (nd == [1*width/8, 0]).all(): dof[i,:] = [1,0];
            if (nd == [7*width/8, 0]).all(): dof[i,:] = [1,0] #Setting the points on the bottom left and bottom right as fixed in y direction
        else:
            if (nd == [0, 0]).all(): dof[i,:] = [0,0];
            if (nd == [width, 0]).all(): dof[i,:] = [1,0] #Setting the points on the bottom left and bottom right as fixed in y direction
        if case == 1 or case == 2 or case == 3: f += [0, -1] if (nd == [width / 2, 0]).all() else [0, 0]
        if case == 4:
            if (nd == [width / 2, height/2]).all(): f += [0, -4]
            elif (nd == [3*width / 8, height/2]).all() or (nd == [5*width / 8, height/2]).all(): f += [0, -1]
            elif (nd == [1*width / 4, height/2]).all() or (nd == [3*width / 4, height/2]).all(): f += [0, -1]
            elif (nd == [1* width / 8, height / 2]).all() or (nd == [7 * width / 8, height / 2]).all(): f += [0, -1]
            elif (nd == [0, height / 2]).all() or (nd == [width, height / 2]).all(): f += [0, -2]
            else: f+= [0,0]
    #Create the 'ground structure'
    for i, j in itertools.combinations(range(len(Nd)), 2): #Sooo cool, getting all the combinations of 2 numbers from the list, aka the i and j node number for all possible elements
        dx, dy = abs(Nd[i][0]-Nd[j][0]), abs(Nd[i][1]-Nd[j][1])
        if gcd(int(dx), int(dy)) == 1 or jc != 0: #Checking if the member would pass through any other points (gcd) or if the joint cost is greater than 0
            seg = [] if convex else LineString([Nd[i], Nd[j]])  #Adding a blank seg vector if the shape is convex, or a 2-point line string with the i and j nodes if shape is not convex
            if convex or poly.contains(seg) or poly.boundary.contains(seg): #Checking if we can add the member
                PML.append([i,j,np.sqrt(dx**2+dy**2), False]) #
    PML, dof = np.array(PML), np.array(dof).flatten() #Setting the member to iterate over in a list PML that contains i node index, y node index, length, and False (no clue why it contains False!)
    f = [f[i:i+len(Nd)*2] for i in range(0, len(f), len(Nd)*2)] #Dont understand this line. For the i value in = [1, ????Still dont get it
    print('Nodes: %d Members: %d' % (len(Nd), len(PML)))
    #Solving
    Cn = PML    # All possible members
    plotTrussGrid(Nd, Cn)
    vol, a, q = solveLP(Nd, Cn, f, dof, st, sc, jc)
    print("Vol: %f, mems: %d" % (vol, len(Cn)))
    print("Volume: %f" % (vol))
    plotTruss(Nd, Cn, a, q, max(a) * 1e-3, (string+"\nVol="+str(round(vol, 5))), False)

##YOU CAN CHANGE THE OPTIONS HERE ##
case = 1 #Change between different cases
text = 1 #Switch to show/hide text on graphs. 1 to show, 0 to hide.

if case == 1:
    trussopt(width=4, height=1, st=1, sc = 1, jc =0, string = 'Case 1 Optimized Truss')
elif case == 2:
    trussopt(width=4, height=1, st=1, sc=1, jc=0, string = 'Case 2 optimized Truss, Max Section Area of .5')
elif case == 3:
    trussopt(width=4, height=2, st=1, sc=1, jc=0, string='Case 3 optimized Truss, Double Height')
elif case == 4:
    trussopt(width=8, height=4, st=1, sc=1, jc=0, string='Case 4 optimized Truss, Double height, No-Go Zone')


