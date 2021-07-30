$$
\omega = \nabla \times u
$$

```
def addVorticity():
    for k in range(1,Nz-1):
        for j in range(1,Ny-1):
            for i in range(1,Nx-1):
                
                right = i + (obstacle[i+1,j,k] == 0)
                left = i - (obstacle[i-1,j,k] == 0)
                up = j + (obstacle[i,j+1,k] == 0)
                down = j - (obstacle[i,j-1,k] == 0)
                front = k + (obstacle[i,j,k+1] == 0)
                back = k - (obstacle[i,j,k-1] == 0)
                
                diffx = max(right - left,1)
                diffy = max(up - down,1)
                diffz = max(front - back,1)
                
                xVort[i,j,k] = (zVel[i,up,k]-zVel[i,down,k])/diffy + (-yVel[i,j,front]+yVel[i,j,back])/diffz
                yVort[i,j,k] = (xVel[i,j,front]-xVel[i,j,back])/diffz + (-zVel[right,j,k]+zVel[left,j,k])/diffx
                zVort[i,j,k] = (yVel[right,j,k]-yVel[left,j,k])/diffx + (-xVel[i,up,k]+xVel[i,down,k])/diffy   
                vorticity[i,j,k] = np.sqrt(xVort[i,j,k]**2 + yVort[i,j,k]**2 + zVort[i,j,k]**2)
```

反过来
$$
\bold u(\bold x_i) = \frac{1}{4\pi} \sum_{j=1,j=i}^n \frac{\omega_j \times (\bold x_i - \bold x_j)}{|\bold x_i - \bold x_j|^3}
$$

```python
		double3 ppos= make_double3(pos[idx].x, pos[idx].y, pos[idx].z);
		if(ppos.x<bbmin.x || ppos.x>bbmax.x 
		|| ppos.y<bbmin.y || ppos.y>bbmax.y
		|| ppos.z<bbmin.z || ppos.z>bbmax.z)
		{
			double3 vel = make_double3(0,0,0);
			double3 dir = make_double3(ppos.x-vort_center.x,ppos.y-vort_center.y,ppos.z-vort_center.z);
			double r3 = compute_r3(dir);
			vel = cross_uxv(omega, dir);
            # 但是这个vel 究竟是网格还是粒子的？
			vel.x = vel.x/r3;
			vel.y = vel.y/r3;
			vel.z = vel.z/r3;

			u[idx]=0;
			v[idx]=0;
			w[idx]=0;
		}
		
def cross_uxv(u,v):
    term0 = u[1]*v[2] - u[2]*v[1] # UyVz - UzVy
    term1 = u[2]*v[0] - u[0]*v[2] # UzVx - UxVz
    term2 = u[0]*v[1] - u[1]*v[0] # UxVy - UyVx
    return np.array([term0,term1,term2])
```

which is known as Biot-Savart summation

Baroclinic Turbulence with Varying Density and Temperature  
$$
\frac{D \omega}{Dt} = (\omega \cdot \nabla)\bold u + \frac{1}{\rho^2}\nabla \rho \times p + \nabla \times \bold B + \nu \nabla ^2 \omega
$$
右边第一项是涡量的strething或者tilting，由速度梯度导致。

右边第二项是斜压项barcolinic。

第三项是体积力，比如重力。Unlike the advection, stretching and tilting terms in equation 8, this baroclinic vector can generate vorticity in cases where the pressure gradient and thermal gradient are not parallel.  