#Minimizing an objective function of the Volumne, Constructibility, and Strain Energy 
 
set MEMBERS; 
set NODES; 
param Tf = 4896000; 
param Smax = 40; #creating a maximum section if we want 
param long {MEMBERS} > 0;

param fx {NODES}; #External force in x 
param fy {NODES}; #External force in y 
param dofx {NODES}; #Whether or not x at the node is a DOF 
param dofy {NODES}; #Whether or not y at the node is a DOF
param Cx  {NODES,MEMBERS}; # 
param Cy {NODES,MEMBERS};
param jc = 1; #Joint cost 
param E = 4176000000; #Young's Modulus 
param Cw = 1; #How much we car about weight
param Cc = 1; #How much we car about constructibility
param Cd = 0; #How much we car about deflection
param Fext = 2600; #Magnitude of the external force

var ID {j in MEMBERS} integer <=1, >=0, :=1;  #Existence of bar
#var ID {j in MEMBERS} binary :=1;  #Existence of bar
var S {j in MEMBERS} >= 0.0000001, <= 40;  # Seccion of each bar
var F {j in MEMBERS};  # Force in each bar
minimize tri_obj: sum {j in MEMBERS} (Cw*(S[j]*(long[j]))+Cc*(ID[j])+(Cc*S[j]*(jc))+Cd*(1/(2*Fext)*((F[j])^2)*long[j]/(S[j]*E)));
subject to Fx {i in NODES}: sum {j in MEMBERS} F[j] * Cx[i,j]*ID[j] = -fx[i]*dofx[i]; #Node force equilibrium in x
subject to Fy {i in NODES}: sum {j in MEMBERS} F[j] * Cy[i,j]*ID[j] = -fy[i]*dofy[i]; #Node force equilibrium in y
subject to max_axial {j in MEMBERS}: F[j]*ID[j] <= Tf*S[j];  #Max tensile stress
subject to min_axial {j in MEMBERS}: F[j]*ID[j] >= -Tf*S[j]; #Max compressive Stress

 data;  ############ DATA STARTS HERE ############
set MEMBERS :=    b1   b2   b3   b4   b5   b6   b7   b8   b9   b10   b11   b12   b13   b14   b15;
set NODES :=  n1 n2 n3 n4 n5 n6;
param:   long :=
 b1   1.0
 b2   2.0
 b3   1.0
 b4   1.4142135623730951
 b5   2.23606797749979
 b6   1.0
 b7   1.4142135623730951
 b8   1.0
 b9   1.4142135623730951
 b10   2.23606797749979
 b11   1.4142135623730951
 b12   1.0
 b13   1.0
 b14   2.0
 b15   1.0;

param : dofx dofy :=
 n1   1   0
 n2   1   1
 n3   1   0
 n4   1   1
 n5   1   1
 n6   1   1;

param : fx fy :=
n1   0   0
n2   0   0
n3   0   0
n4   0   0
n5   0   -2600
n6   0   0;

param Cx :
       b1   b2   b3   b4   b5   b6   b7   b8   b9   b10   b11   b12   b13   b14   b15 :=
n1   -1.0   -1.0   0.0   -0.7071067811865475   -0.8944271909999159   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n2   1.0   0.0   0.0   0.0   0.0   -1.0   0.7071067811865475   0.0   -0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0
n3   0.0   1.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.8944271909999159   0.7071067811865475   0.0   0.0   0.0   0.0
n4   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   0.0   0.0   -0.8944271909999159   0.0   0.0   -1.0   -1.0   0.0
n5   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   0.0   1.0   0.0   -1.0
n6   0.0   0.0   0.0   0.0   0.8944271909999159   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   1.0   1.0;

param Cy :
       b1   b2   b3   b4   b5   b6   b7   b8   b9   b10   b11   b12   b13   b14   b15 :=
n1   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n2   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   -1.0   -0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0
n3   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n4   0.0   0.0   1.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.0   0.0
n5   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   1.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0
n6   0.0   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   1.0   0.0   0.0   0.0;