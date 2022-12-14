#Minimizing an objective function of the number of members  ##Volumne, Constructibility, and Strain Energy 
 
set MEMBERS; 
set NODES; 
param Tf = 4896000; 
param Smax = 0.5; #creating a maximum section if we want 
param long {MEMBERS} > 0;

param fx {NODES}; #External force in x 
param fy {NODES}; #External force in y 
param dofx {NODES}; #Whether or not x at the node is a DOF 
param dofy {NODES}; #Whether or not y at the node is a DOF
param Cx  {NODES,MEMBERS}; # 
param Cy {NODES,MEMBERS};
param jc = 0; #Joint cost 
param E = 4176000000; #Young's Modulus 
param Cw = 0; #How much we car about weight
param Cc = 1000.0; #How much we car about constructibility
#param Cd = 0; #How much we car about deflection
param Fext = 1300; #Magnitude of the external force

var ID {j in MEMBERS} integer >=0, <=1, :=1;  #Existence of bar
var S {j in MEMBERS} >= 0.00000001, <= 0.5, :=.1;  # Seccion of each bar
var F {j in MEMBERS}, :=2600;  # Force in each bar
minimize tri_obj: sum {j in MEMBERS} (.0001*Cc*(S[j]*(long[j]))+Cc*(ID[j]));
subject to Fx {i in NODES}: sum {j in MEMBERS} F[j] * Cx[i,j]*ID[j] = -fx[i]*dofx[i]; #Node force equilibrium in x
subject to Fy {i in NODES}: sum {j in MEMBERS} F[j] * Cy[i,j]*ID[j] = -fy[i]*dofy[i]; #Node force equilibrium in y
subject to max_axial {j in MEMBERS}: F[j]*ID[j] <= Tf*S[j];  #Max tensile stress
subject to min_axial {j in MEMBERS}: F[j]*ID[j] >= -Tf*S[j]; #Max compressive Stress

 data;  ############ DATA STARTS HERE ############
set MEMBERS :=    b1   b2   b3   b4   b5   b6   b7   b8   b9   b10   b11   b12   b13   b14   b15   b16   b17   b18   b19   b20   b21   b22   b23   b24   b25   b26   b27   b28   b29   b30   b31   b32   b33   b34   b35   b36   b37   b38   b39   b40   b41;
set NODES :=  n1 n2 n3 n4 n5 n6 n7 n8 n9 n10;
param:   long :=
 b1   1.0
 b2   2.0
 b3   3.0
 b4   1.0
 b5   1.4142135623730951
 b6   2.23606797749979
 b7   3.1622776601683795
 b8   1.0
 b9   2.0
 b10   3.0
 b11   1.4142135623730951
 b12   1.0
 b13   1.4142135623730951
 b14   2.23606797749979
 b15   3.1622776601683795
 b16   1.0
 b17   2.0
 b18   2.23606797749979
 b19   1.4142135623730951
 b20   1.0
 b21   1.4142135623730951
 b22   2.23606797749979
 b23   1.0
 b24   3.1622776601683795
 b25   2.23606797749979
 b26   1.4142135623730951
 b27   1.0
 b28   1.4142135623730951
 b29   3.1622776601683795
 b30   2.23606797749979
 b31   1.4142135623730951
 b32   1.0
 b33   1.0
 b34   2.0
 b35   3.0
 b36   1.0
 b37   2.0
 b38   3.0
 b39   1.0
 b40   2.0
 b41   1.0;

param : dofx dofy :=
 n1   1   0
 n2   1   1
 n3   1   1
 n4   1   1
 n5   1   0
 n6   1   1
 n7   1   1
 n8   1   1
 n9   1   1
 n10   1   1;

param : fx fy :=
n1   0   0
n2   0   0
n3   0   -1300
n4   0   0
n5   0   0
n6   0   0
n7   0   0
n8   0   0
n9   0   0
n10   0   0;

param Cx :
       b1   b2   b3   b4   b5   b6   b7   b8   b9   b10   b11   b12   b13   b14   b15   b16   b17   b18   b19   b20   b21   b22   b23   b24   b25   b26   b27   b28   b29   b30   b31   b32   b33   b34   b35   b36   b37   b38   b39   b40   b41 :=
n1   -1.0   -1.0   -1.0   0.0   -0.7071067811865475   -0.8944271909999159   -0.9486832980505138   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n2   1.0   0.0   0.0   0.0   0.0   0.0   0.0   -1.0   -1.0   -1.0   0.7071067811865475   0.0   -0.7071067811865475   -0.8944271909999159   -0.9486832980505138   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n3   0.0   1.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -1.0   -1.0   0.8944271909999159   0.7071067811865475   0.0   -0.7071067811865475   -0.8944271909999159   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n4   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   -1.0   0.9486832980505138   0.8944271909999159   0.7071067811865475   0.0   -0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n5   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.9486832980505138   0.8944271909999159   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n6   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   -0.8944271909999159   0.0   0.0   0.0   0.0   0.0   -0.9486832980505138   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -1.0   -1.0   -1.0   0.0   0.0   0.0   0.0   0.0   0.0
n7   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   0.0   0.0   0.0   0.0   0.0   -0.8944271909999159   0.0   0.0   0.0   -0.9486832980505138   0.0   0.0   0.0   1.0   0.0   0.0   -1.0   -1.0   -1.0   0.0   0.0   0.0
n8   0.0   0.0   0.0   0.0   0.0   0.8944271909999159   0.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   0.0   0.0   0.0   -0.8944271909999159   0.0   0.0   0.0   1.0   0.0   1.0   0.0   0.0   -1.0   -1.0   0.0
n9   0.0   0.0   0.0   0.0   0.0   0.0   0.9486832980505138   0.0   0.0   0.0   0.0   0.0   0.0   0.8944271909999159   0.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   0.0   0.0   0.0   1.0   0.0   1.0   0.0   1.0   0.0   -1.0
n10   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.9486832980505138   0.0   0.0   0.0   0.0   0.0   0.0   0.8944271909999159   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   1.0   1.0;

param Cy :
       b1   b2   b3   b4   b5   b6   b7   b8   b9   b10   b11   b12   b13   b14   b15   b16   b17   b18   b19   b20   b21   b22   b23   b24   b25   b26   b27   b28   b29   b30   b31   b32   b33   b34   b35   b36   b37   b38   b39   b40   b41 :=
n1   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n2   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.7071067811865475   -1.0   -0.7071067811865475   -0.4472135954999579   -0.31622776601683794   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n3   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.4472135954999579   -0.7071067811865475   -1.0   -0.7071067811865475   -0.4472135954999579   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n4   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   -0.31622776601683794   -0.4472135954999579   -0.7071067811865475   -1.0   -0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n5   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n6   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.0   0.0   0.31622776601683794   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n7   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.31622776601683794   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n8   0.0   0.0   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n9   0.0   0.0   0.0   0.0   0.0   0.0   0.31622776601683794   0.0   0.0   0.0   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
n10   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.31622776601683794   0.0   0.0   0.0   0.0   0.0   0.0   0.4472135954999579   0.0   0.0   0.0   0.0   0.0   0.7071067811865475   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0;