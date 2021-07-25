import numpy as np

q = 1
k = 1
nseg = 8
rq = 1


Npoints = 39

coord = np.array([[0,1.1],[0,1.5],[0,2],[0,2.5],[0,3],[0,3.5],[0,4],[0,4.5],
                  [0,5],[0,-1.1],[0,-1.5],[0,-2],[0,-3],[0,-3.5],[0,-4],
                  [0,-4.5],[0,-5],[1.1,0],[1.15,0],[1.3,0],[1.5,0],
                  [2,0],[2.5,0],[3,0],[3.5,0],[4,0],[4.5,0],[5,0],
                  [5,5,0],[6,0],[6.5,0],[7,0],[7.5,0],[8,0],[8.5,0],
                  [9,0],[9.5,0],[10,0]])

C = 0.5 / np.pi
cl = np.pi * k / 2
delta = 2 * np.pi / nseg # increment in angle theta
rqo = rq / np.cos(delta / 2) # outer radius of isolator
theta = (np.pi - delta) / 2

xA = np.zeros((nseg,2))
xB = np.zeros((nseg,2))
xS = np.zeros((nseg,2))

xA[0,0] = rqo * np.cos(theta)
xA[0,1] = rqo * np.sin(theta)

for i in range(nseg-1):
    xB[i,0] = rqo * np.cos(theta)
    xB[i,1] = rqo * np.sin(theta)
    xA[i+1,0] = xB[i,0]
    xA[i+1,1] = xB[i,1]
    
xB[nseg-1,0] = xA[0,0]
xB[nseg-1,1] = xA[0,1]

for i in range(nseg):
    xS[i,0] = (xA[i,0] + xB[i,0]) / 2
    xS[i,1] = (xA[i,1] + xA[i,1]) / 2
    
theta = np.pi / 2
t0 = np.zeros((nseg))
for i in range(nseg):
    t0[i] = q * np.sin(theta)
    theta = theta + delta
    
# Assemble matrices DT and DU
ve = np.zeros((nseg,2))
vn = np.zeros((nseg,2))
Lhs = np.zeros((nseg,nseg))
Rhs = np.zeros((nseg,nseg))

for i in range(nseg):
    dx = xA[i,0] - xB[i,0]
    dy = xA[i,1] - xB[i,1]
    dis = np.sqrt(dx*dx + dy*dy)
    ve[i,0] = dx / dis
    ve[i,1] = dy / dis
    vn[i,0] = ve[i,1]
    vn[i,1] = -ve[i,0]
    for j in range(nseg):
        rA = np.sqrt((xA[i,0]-xS[i,0])**2 + (xA[i,1]-xS[i,1])**2)
        rB = np.sqrt((xB[i,0]-xS[i,0])**2 + (xB[i,1]-xS[i,1])**2)
        vrA = [xA[nseg,0] - xS[nseg,0],xA[nseg,1] - xS[nseg,1]]
        vrB = [xB[nseg,0] - xS[nseg,0],xB[nseg,1] - xS[nseg,1]]
        cosThetaA = - np.dot(vn[nseg,:],vrA) / rA
        cosThetaB = - np.dot(vn[nseg,:],vrB) / rB
        sinThetaA = - np.dot(ve[nseg,:],vrA) / rA
        sinThetaB = - np.dot(ve[nseg,:],vrB) / rB
        thetaA = np.acos(cosThetaA) * np.sign(sinThetaA)
        thetaB = np.acos(cosThetaB) * np.sign(sinThetaB)
        if i == j:
            Lhs[i,i] = 5
            Rhs[i,i] = dis * cl * (np.log(dis / 2) - 1)
        else:
            Lhs[i,j] = C * (thetaB - thetaA)
            Rhs[i,j] = cl*(rB*sinThetaB*(np.log(rB)-1)+thetaB*rB*cosThetaB
                        -rA*sinThetaA*(np.log(rA)-1)-thetaA*rA*cosThetaA)
            
F = np.dot(Rhs,t0)
u = np.dot(np.linalg.inv(Lhs),F)
for i in range(Npoints):
    up = 0
    qx = 0
    qy = 0
    for j in range(nseg):
        rA = np.sqrt((xA[i,0]-xS[i,0])**2 + (xA[i,1]-xS[i,1])**2)
        rB = np.sqrt((xB[i,0]-xS[i,0])**2 + (xB[i,1]-xS[i,1])**2)
        vrA = [xA[nseg,0] - xS[nseg,0],xA[nseg,1] - xS[nseg,1]]
        vrB = [xB[nseg,0] - xS[nseg,0],xB[nseg,1] - xS[nseg,1]]
        cosThetaA = - np.dot(vn[nseg,:],vrA) / rA
        cosThetaB = - np.dot(vn[nseg,:],vrB) / rB
        sinThetaA = - np.dot(ve[nseg,:],vrA) / rA
        sinThetaB = - np.dot(ve[nseg,:],vrB) / rB
        H = rA * cosThetaA
        thetaA = np.acos(cosThetaA) * np.sign(sinThetaA)
        thetaB = np.acos(cosThetaB) * np.sign(sinThetaB)
        if thetaB - thetaA > np.pi:
            thetaA = 2 * np.pi
        dT = C * (thetaB - thetaA)
        dU = cl * (rB*sinThetaB*(np.log(rB)-1)+thetaB*rB*cosThetaB
                        -rA*sinThetaA*(np.log(rA)-1)-thetaA*rA*cosThetaA)
        dSx = C / k * (thetaB - thetaA)
        Fact = cosThetaB / cosThetaA
        dSy = 0
        if Fact > 0:
            dSy = - C / k * np.log(Fact)
        dRx = - C / H * (cosThetaB*sinThetaB - cosThetaA*sinThetaA)
        dRy = C / H * (cosThetaB**2 - cosThetaA**2)
        up = up + dU*t0[j] - dT*u[j]
        qxp = -k * (dSx*t0[j] - dRx*u[j])
        qyp = -k * (dSy*t0[j] - dRy*u[j])
        qx = qx + qxp*vn[j,0] - qyp*vn[j,1]
        qy = qy + qxp*vn[j,1] + qyp*vn[j,0]
    Up = Up - q / k * xS[i,1]
    qy = qy + q
        