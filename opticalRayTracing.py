# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.axes as axs
import scipy
import math
# %matplotlib inline

def RayMatrix(n,angle,height):
    return np.array(((n*angle), (height)))

def RefractionMatrix(n_incident,n_transmitted,R):
    return np.array([[1,(-1*(n_transmitted-n_incident))/R], [0,1]])

def TransferMatrix(n,d):
    return np.array([[1,0],[d/n , 1]])

#CONSTANTS==================================================================================
    d_object = 20.0                 #distance of object from first lens 
    y_object = [-10,10]             #height of object

    r1 = 5.0                        #radius of curvature of first lens, initial side
    r2 = -5.0                       #radius of curvature of first lens, transmitted side
    r3 = 10.0                       #radius of curvature of second lens, initial side
    r4 = -10.0                      #radius of curvature of second lens, transmitted side

    n1 = 1.0                        #index of refraction of initial medium
    n2 = 1.5                        #index of refraction of first lens
    n3 = 1.0                        #index of refraction after the first lens
    n4 = 1.5                        #index of refraction of second lens
    n5 = 1.0                        #index of refraction after the second lens

    d1 = 2.0                        #thickness of first lens
    d2 = 2.0                        #thickness of second lens
    
    l = 12                          #distance between two lens
    l2 = 40                         #arbitrary distance of open air after second lens

def matrix_system(height):
  yi_1 = height  
  angle = 0

  r1_incident = RayMatrix(n1,angle,yi_1)
  R1 = RefractionMatrix(n1,n2,r1)
  r1_transmitted = R1 @ r1_incident 

  T21 = TransferMatrix(n2,d1)
  r2_incident = T21 @ r1_transmitted

  R2 = RefractionMatrix(n2,n3,r2)
  r2_transmitted = R2 @ r2_incident

  A1 = R2 @ (T21 @ R1)
  
  T_open = TransferMatrix(n3,l)
  r3_incident = T_open @ r2_transmitted

  R3 = RefractionMatrix(n3,n4,r3)
  r3_transmitted = R3 @ r3_incident

  T43 = TransferMatrix(n4,d2)
  r4_incident = T43 @ r3_transmitted

  R4 = RefractionMatrix(n4,n5,r4)
  r4_transmitted = R4 @ r4_incident
  
  A2 = R4 @ (T43 @ R3)

  T56 = TransferMatrix(n5,l2)
  r5_incident = T56 @ r4_transmitted
  
  return r1_incident,R1,r1_transmitted,T21,r2_incident,R2,r2_transmitted,T_open,r3_incident,R3,r3_transmitted,T43,r4_incident,R4,r4_transmitted,T56,r5_incident,A1,A2

def principal_planes(A,d1,d2):
    #FOCAL LENGTH AND PRINCIPAL PLANES===========================================================
    fo = 1/A[0][1]                         #focal length on the initial side
    fi = -1*fo                             #focal length on the transmitted side
    v1h1 = (n1*(1 - A[0][0]))/(-1*A[0][1]) #first principal plane of lens
    v2h2 = (n3*(A[1][1] - 1))/(-1*A[0][1]) #second principal plane of lens
    
    #DETERMINING FRONT FOCAL AND BACK FOCAL========================================================
    frontfocal = A[0][0]*fo                                     #front focus
    backfocal = A[1][1]*fi                                      #back focus
    #effec_focal = (n2-1)*((1/r1)-(1/r2)+((n2-1)*d)/(r1*r2*n2)) #effective focal length

    #PLOT PRINCIPAL PLANES AND FOCII========================================================
    plt.axvline(d1+v1h1, color='red')
    plt.axvline(d2+v2h2, color='red')
    plt.plot(d1+fo,0,'bo')
    plt.plot(d2+fi,0, 'go')

def image_formed(A,yi_1,d_object):
    angle =  np.arctan(yi_1/d_object)
    d_Image = (-1*(A[1][0]) - (A[1][1])*d_object)/((A[0][0])+(A[0][1])*d_object)                                             #distance of the image from the lens
    y_Image =angle*(A[1][0] + A[1][1]*d_object + (A[0][0]+A[0][1]*d_object)*d_Image)+ yi_1*(A[1][1] + A[0][1]*d_Image) #height of image from the optical plane
    return d_Image,y_Image

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

plt.figure(figsize =(20,7) )

#PLOTTING THE OBJECT=================================================
plt.plot([-1*d_object,-1*d_object],[min(y_object),max(y_object)], color='red')

#PLOTTING THE PRINCINPAL PLANES =====================================
plt.axvline(0, color='black')    #first surface of first lens
plt.axvline(d1, color='black')    #second surface of first lens
plt.axvline(d1+l, color='black')    #first surface of second lens
plt.axvline(d1+l+d2, color='black')    #second surface of second lens
plt.axhline(0, color='black')    #optical axis

y1,y2 = [], []
#RAY TRACING =================================================================================
for r in y_object:
    r1_incident,R1,r1_transmitted,T21,r2_incident,R2,r2_transmitted,T_open,r3_incident,R3,r3_transmitted,T43,r4_incident,R4,r4_transmitted,T56,r5_incident,A1,A2 = matrix_system(r) 

    principal_planes(A1,0,d1)    
    
    plt.plot([-1*d_object,0],[r,r1_incident[1]], 'b:') #ray from object to P1
    plt.plot([0,d1],[r1_transmitted[1],r2_incident[1]], 'b:')  #refracted ray through the first half of the first lens
    plt.plot([d1,d1+l],[r2_transmitted[1],r3_incident[1]], 'b:') #ray from first lens to second lens
    
    principal_planes(A2,d1+l,d1+l+d2)

    plt.plot([d1+l,d1+l+d2],[r3_transmitted[1],r4_incident[1]], 'b:') #refracted ray through the first half of the second lens
    plt.plot([d1+l+d2,d1+l+d2+l2],[r4_transmitted[1],r5_incident[1]], 'b:') #ray from second lens onwards
    y1.append(r4_transmitted[1])
    y2.append(r5_incident[1])
    

x = [d1+l+d2,d1+l+d2+l2]

L1 = line([x[0],y1[0]], [x[1],y2[0]])
L2 = line([x[0],y1[1]], [x[1],y2[1]])
R = intersection(L1,L2)

if R != False:
    plt.plot([x[0],R[0]],[y1[0],R[1]],'g:')
    plt.plot([x[0],R[0]],[y1[1],R[1]],'g:')
    plt.plot(R[0],R[1], 'ro')
plt.savefig('case1.png')
plt.show()


