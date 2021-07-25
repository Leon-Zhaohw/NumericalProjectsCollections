import numpy as np
import Laplace
'''
The Boundary Element 
Trefftz Method

仅供代码参考，结果可能不正确

'''

q = 1 # Heat inflow
k = 1 # Thermal conductivity
npnts = 32 # Number of points P,Q
rq = 1 # radius of isolator
rp = 0.7 # radius of sources

Npoints = 18

coord = np.array([[0,-5],[0,-4.5],[0,-4],
                  [0,-3.5],[0,-3],[0,-2.5],
                  [0,-2],[0,-1.5],[0,-1],
                  [0,1],[0,1.5],[0,2],
                  [0,2.5],[0,3],[0,3.5],
                  [0,4],[0,4.5],[0,5]])

Delta = 2 * np.pi / npnts # increment in angle theta between points
thetaQ = np.pi / 2 # angle theta to first field point Q1

Rhs = np.zeros((npnts))
Lhs = np.zeros((npnts,npnts))


for i in range(npnts): # Field Point
    Rhs[i] = q * np.sin(thetaQ)
    xq = rq * np.cos(thetaQ)
    yq = rq * np.sin(thetaQ)
    
    vnorm = [-np.cos(thetaQ),-np.sin(thetaQ)]
    
    thetaQ += Delta
    thetaP = np.pi / 2
    
    # 这么多个在j处的Source point，对在 i 处的Field Point，
    # 距离为 r，造成的影响为的比例为Lhs[i,j]
    # 在 j 处的场点，它的强度为F[j]，它的距离带来的影响为Lhs[:,j]
    
    for j in range(npnts): # Source Point
        xp = rp * np.cos(thetaP)
        yp = rp * np.sin(thetaP)
        dxr = [xp-xq,yp-yq]
        r = np.sqrt(dxr[0]**2 + dxr[1]**2)
        dxr /= r
        Lhs[i,j] = Laplace.T(r,dxr,vnorm,2)
        thetaP += Delta
    
Lhs = - Lhs
Rhs = - Rhs  
F = np.dot(np.linalg.inv(Lhs),Rhs)
thetaQ = np.pi / 2
  
# 问题1，uq和Rhs的区别是什么，为什么U和T不一样

for i in range(npnts): # Field Point
    uq = 0
    xq = rq * np.cos(thetaQ)
    yq = rq * np.sin(thetaQ)
    
    thetaQ = thetaQ + Delta
    thetaP = np.pi / 2
    
    for j in range(npnts): # Source Point
        xp = rp * np.cos(thetaP)
        yp = rp * np.sin(thetaP)
        
        dxr = [xp-xq,yp-yq]
        
        r = np.sqrt(dxr[0]**2 + dxr[1]**2)
        uq = uq + Laplace.U(r, k,2)*F[j]
        thetaP += Delta
        
    uq = uq - q / k * yq
    print("Temprature at field point ",i," = ",uq)
    # temperature at field point uq
    
for i in range(Npoints):
    uq = 0
    qx = 0
    qy = 0
    thetaP = np.pi / 2
    xi = coord[i,0]
    yi = coord[i,1]
    for j in range(npnts):
        xp = rp * np.cos(thetaP)
        yp = rp * np.sin(thetaP)
        dx = [-xp + xi, -yp + yi]
        r = np.sqrt(dx[0]**2 + dx[1]**2)
        dxr = dx / r
        uq = uq + Laplace.U(r,k,2)*F[j]
        dUxy = Laplace.dU(r,dxr,2)
        qx = qx + dUxy[0]*F[j]
        qy = qy + dUxy[1]*F[j]
        thetaP = thetaP + Delta
    uq = uq - q / k * yi
    qy = qy + 1
    print("at x = ",xi,"y = ",yi,"temprature = ",uq," qx = ",qx,"qy = ",qy)
