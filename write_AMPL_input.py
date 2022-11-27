from math import gcd, ceil
import itertools
from scipy import sparse
import pandas as pd
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

def trussopt_settup(width, height, st, sc, jc):
    #Creating the original polygon shape
    case = 1
    if case == 1 or case == 2:
        poly = Polygon([(0,0), (width, 0), (width, height), (0, height)]) #Create the area over which to iterate
    elif case == 3:
        ext = [(0, 0), (width, 0), (width, height), (0, height), (0, 0)]
        interior = [(width/2-5, 1), (width/2+5, 1), (width/2+5, 7), (width/2-5, 7), (width/2-5, 1)][::-1]
        poly = Polygon(ext, [interior])
    convex = True if poly.convex_hull.area == poly.area else False  #checks weather area is convex
    xv, yv = np.meshgrid(range(width+1), range(height+1))  #defines a meshgrid with the possible points
    pts = [Point(xv.flat[i], yv.flat[i]) for i in range(xv.size)] #crete a list of points for the possible nodes
    Nd = np.array([[pt.x,pt.y] for pt in pts if poly.intersects(pt)]) #Creating a list of the possible nodes (note the the "if" statement excludes points that would be within a central nogo zone
    dof, f, PML = np.ones((len(Nd), 2)), [], [] #Setting a dof vector to be an array of ones, NX2, also defining f vec and PML vec (not sure what these are?)
    #Load and support conditions
    for i, nd in enumerate(Nd): #Looping through all the nodes one by one
        if (nd == [0, 0]).all() or (nd == [width, 0]).all() : dof[i,:] = [1,0] #Setting the points on the bottom left and bottom right as fixed in y direction
        # #Case 1
        if case == 1 or case == 2: f += [0, -1] if (nd == [width / 2, 0]).all() else [0, 0]
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
    Cn = PML #Only iterating over members with lengths shorter or equal to root 2
    B = calcB(Nd, Cn, dof)
    plotTrussGrid(Nd, Cn)
    return Nd, Cn, f, dof, B


def write2Ampl(Tf,Smax,jc,E,Cw,Cc,Cd,Cn,Fext,f,dof,B,ampl_model,ampl_run,ampl_NEOS_run):
    with open(ampl_model, "w") as inp:
        #Intro and set parameters
        inp.write("#Minimizing an objective function of the Volumne, Constructibility, and Strain Energy \n \nset MEMBERS; \nset NODES; \nparam Tf = "+str(1)+"; \nparam Smax = "+str(Smax)+"; #creating a maximum section if we want \nparam long {MEMBERS} > 0;")
        inp.write("\n\nparam fx {NODES}; #External force in x \nparam fy {NODES}; #External force in y ")
        inp.write("\nparam dofx {NODES}; #Whether or not x at the node is a DOF \nparam dofy {NODES}; #Whether or not y at the node is a DOF")
        inp.write("\nparam Cx  {NODES,MEMBERS}; # \nparam Cy {NODES,MEMBERS};")
        inp.write("\nparam jc = "+str(jc)+"; #Joint cost \nparam E = "+str(E)+"; #Young's Modulus ")
        inp.write("\nparam Cw = "+str(Cw)+"; #How much we car about weight")
        inp.write("\nparam Cc = "+str(Cc)+"; #How much we car about constructibility")
        inp.write("\nparam Cd = " + str(Cd) + "; #How much we car about deflection")
        inp.write("\nparam Fext = " + str(Fext) + "; #Magnitude of the external force")
        #Defining variables
        inp.write("\n\nvar S {j in MEMBERS} >= 0.00001, <= "+str(Smax)+";  # Seccion of each bar")
        inp.write("\nvar F {j in MEMBERS};  # Force in each bar")
        #Ojectve function
        # inp.write("\nminimize tri_obj: sum {j in MEMBERS} (S[j]*(long[j]+jc)+Cd*((F[j])^2)*long[j]/(S[j]*E));")
        inp.write("\nminimize tri_obj: sum {j in MEMBERS} (Cw*S[j]*(long[j])+(Cc*S[j]*(jc))+Cd*1/(2*Fext)*((F[j])^2)*long[j]/(S[j]*E));")

        #Constraints
        inp.write("\nsubject to Fx {i in NODES}: sum {j in MEMBERS} F[j] * Cx[i,j] = -fx[i]*dofx[i]; #Node force equilibrium in x")
        inp.write("\nsubject to Fy {i in NODES}: sum {j in MEMBERS} F[j] * Cy[i,j] = -fy[i]*dofy[i]; #Node force equilibrium in y")
        inp.write("\nsubject to max_axial {j in MEMBERS}: F[j] <= Tf*S[j];  #Max tensile stress")
        inp.write("\nsubject to min_axial {j in MEMBERS}: F[j] >= -Tf*S[j]; #Max compressive Stress")
        #DATA
        inp.write("\n\n data;  ############ DATA STARTS HERE ############")
        #Member list
        memlist = []
        memstring = ""
        Nm = int(np.size(Cn)/4)
        for b in range(Nm):
            memlist.append("b"+str(b+1))
            memstring = memstring+"   b"+str(b+1)
        inp.write("\nset MEMBERS := "+memstring+";")
        #Node list
        N = int(np.size(f)/2)
        nodelist = []
        nodestring = ""
        for n in range(N):
            nodelist.append("n"+str(n+1))
            nodestring = nodestring+" n"+str(n+1)
        inp.write("\nset NODES := " + nodestring+";")
        #Longitud de los miembros
        memlenstring = ""
        for b in range(Nm):
            memlenstring = memlenstring+"\n b"+str(b+1)+"   "+str(Cn[b][2])
        inp.write("\nparam:   long :=")
        inp.write(memlenstring+";")
        #dof in x and y direction
        inp.write("\n\nparam : dofx dofy :=")
        dofstring = ""
        for n in range(N):
            dofstring = dofstring+"\n n"+str(n+1)+"   "+str(int(dof[2*n]))+"   "+str(int(dof[2*n+1]))
        inp.write(dofstring+";")
        #External forces
        inp.write("\n\nparam : fx fy :=")
        extforcestring = ""
        for n in range(N):
            extforcestring = extforcestring+"\nn"+str(n+1)+"   "+str(f[0][2*n])+"   "+str(f[0][2*n+1])
        inp.write(extforcestring+";")
        #Cx
        inp.write("\n\nparam Cx :\n    "+memstring+" :=")
        Cx_str = ""
        B_np = B.toarray()
        for n in range(N):
            Cx_str = Cx_str+"\nn"+str(n+1)
            for b in range(Nm):
                Cx_str = Cx_str+"   "+str(B_np[2*n][b])
        inp.write(Cx_str+";")
        #Cy
        inp.write("\n\nparam Cy :\n    " + memstring + " :=")
        Cy_str = ""
        for n in range(N):
            Cy_str = Cy_str + "\nn" + str(n + 1)
            for b in range(Nm):
                Cy_str = Cy_str + "   " + str(B_np[2 * n+1][b])
        inp.write(Cy_str+";")

    with open(ampl_run, "w") as run:
        run.write("reset; model input13.mod; solve; print solve_message >> obj.out; close obj.out; csvdisplay F,S >Trial.csv; close Trial.csv;")

    with open(ampl_NEOS_run, "w") as run:
        run.write(" option solver minos;")
        run.write("\n option show_stats 1;")

        run.write("\n\noption minos_options  ' \ ")
        run.write("\n   superbasics_limit=3000\ ")
        run.write("\n';")
        run.write("\n\nsolve;  \ndisplay F,S;  ")

def AMPL_files_write(jc, Cd, text, case, Tf, Fext, Smax, E, width, height, preprocess, postprocess, ampl_model, ampl_run,ampl_NEOS_run):
    if preprocess == 1:
        Nd, Cn, f, dof, B = trussopt_settup(width=width, height=height, st=1, sc=1, jc=1)
        write2Ampl(Tf, Smax, jc, E, Cw, Cc, Cd, Cn, Fext, f, dof, B, ampl_model, ampl_run, ampl_NEOS_run)

    if postprocess == 1:
        memN = len(Cn)
        df = pd.read_csv(r"C:\Users\saulc\OneDrive\Documents\Optimization\InputTrials\Trial.csv")
        members, area, force = df['Index_1'], df['S'], df['F']
        members_from, map, a, q = np.zeros((memN), dtype=int), np.zeros((memN), dtype=int), np.zeros((memN)), np.zeros(
            (1, memN))
        for i in range(memN):
            members_from[i] = members.iloc[i][1:]
        for i in range(memN):
            a[members_from[i] - 1] = area[i];
            q[0, members_from[i] - 1] = force[i]
        plotTruss(Nd, Cn, a, q, max(a) * 1e-3, string="Optimized Truss", update=True)

jc = 1
Cw = 1
Cc = 1
Cd =100
text = 1  # Switch to show/hide text on graphs. 1 to show, 0 to hide.
case = 1
Tf = 1
Smax = 10
E = 10
width = 20
height = 2
Fext = 1

# to perform post-processing after obtaining an AMPL output code
preprocess = 1
postprocess = 1

#File names and locations
ampl_model = r"C:\Users\saulc\OneDrive\Documents\Optimization\InputTrials\input13.mod"
ampl_run = r"C:\Users\saulc\OneDrive\Documents\Optimization\InputTrials\TrussOptiAMPL.run"
ampl_NEOS_run = r"C:\Users\saulc\OneDrive\Documents\Optimization\InputTrials\NEOSTrussOptiAMPL.run"

#Run main function to create the .mod file and the two .run files
AMPL_files_write(jc, Cd, text, case, Tf, Fext, Smax, E, width, height, preprocess, postprocess, ampl_model, ampl_run,ampl_NEOS_run)

