/* 
Tyler de Witt
2012-08-03
*/

import java.lang.Math;

public class LE {
    // Velocity grid size
    public int mx;
    public int my;

    // Velocity basis fields
    public double [][][][]  vel_basis;

    // Velocity basis coefficients
    public double [] coef;

    // Current velocity field;
    public double [][][] vfield;

    // Basis field eigenvalues
    public double [] eigs;
    public double [] eigs_inv;
    public double [] eigs_inv_root;

    // Structure coefficient matrices
    public SparseMatrix []Ck;
    
    // Viscosity
    public double visc = 0.0;

    // Timestep
    public double dt=0.1;

    // Timestep adjustment factor for particle advection
    public double pdt_mult = 1.0;

    // Small numerical offset for keeping particles within bounds
    public double margin = 1e-7;

    // Basis dimension
    public int N;
    public int N_sqrt;

    // Tables for basis field lookup from vector eigenvalue k1,k2
    public int [][]basis_lookup_table;
    public int [][]basis_rlookup_table;

    // Particles
    public int [] particle_index;
    public double [][][] particles;
    public int num_particles=0;
    public int particle_tail_length=3;

    // Density field

    public double [][]density;

    // Density field size
    public int dmx;
    public int dmy;
    
    // External forces
    public double[] forces_dw;
    public boolean forces_pending;

    LE() {
	this(32,16,32);
    }

    LE(int grid_res, int N,int density_grid_res) {
	this.mx = grid_res;
	this.my = grid_res;
	this.dmx = density_grid_res;
	this.dmy = density_grid_res;

	this.N = N;

	this.vfield = new double[2][this.mx+1][this.my+1];
	this.density = new double[this.dmx][this.dmy];

	this.coef = new double[N];
	this.forces_dw = new double[N];

	System.out.println("Filling lookup table.");
	this.fill_lookup_table();

	System.out.println("Precomputing Basis Fields.");
	this.precompute_basis_fields();

	System.out.println("Precomputing Dynamics.");
	this.precompute_dynamics();

	System.out.println("Done.");
    }


    public void step() {
	//Advance the simulation
	double [] dw = new double[this.N];

	// Calculate current energy
	double prev_e = cur_energy();


	// First order Explicit Euler */
	/*
	for (int k=0;k<this.N;k++)
	    // Calculate C_k matrix vector products
	    dw[k] = this.dot(coef,this.Ck[k].mult(coef));
	*/

	/* RK4: higher order explicit integrator */
	double [][] dwt = new double[4][this.N];
	double [][] qn = new double[4][this.N];

	qn[0] = this.coef;

	for (int k=0;k<this.N;k++) {
	    // Calculate C_k matrix vector products
	    dwt[0][k] = this.dot(qn[0],this.Ck[k].mult(qn[0]));
	    qn[1][k] = qn[0][k] + 0.5*this.dt * dwt[0][k];
	}

	for (int k=0;k<this.N;k++) {
	    dwt[1][k] = this.dot(qn[1],this.Ck[k].mult(qn[1]));
	    qn[2][k] = qn[0][k] + 0.5*this.dt * dwt[1][k];	
	}

	for (int k=0;k<this.N;k++) {
	    dwt[2][k] = this.dot(qn[2],this.Ck[k].mult(qn[2]));
	    qn[3][k] = qn[0][k] + this.dt * dwt[2][k];	
	}

	for (int k=0;k<this.N;k++) {
	    dwt[3][k] = this.dot(qn[3],this.Ck[k].mult(qn[3]));
	    dw[k] = (dwt[0][k] + 2*dwt[1][k] + 2*dwt[2][k] + dwt[3][k])/6.0;
	}	
	
	// Take the explicit step
	for (int k=0;k<this.N; k++) 
	    this.coef[k] += dw[k] * this.dt;

	// Renormalize energy
	if (prev_e > 1e-5)
	    this.set_energy(prev_e);


	// Dissipate energy for viscosity
	for (int k=0;k<this.N;k++)  {
	    coef[k] *= Math.exp(-1.0 * this.eigs[k] * this.dt * this.visc);
	    // Add external forces 
	    coef[k] += this.forces_dw[k];
	    forces_dw[k] = 0;
	}

	// Reconstruct velocity field
	this.expand_basis();
    }

    public void precompute_basis_fields() {
	this.vel_basis = new double[this.N][][][];

	for (int i=0;i<this.N;i++) {
	    int k1=this.basis_lookup(i,0);
	    int k2=this.basis_lookup(i,1);

	    this.vel_basis[i] = this.basis_field_2d_rect(k1,k2,1.0);
	}
    }

    public void precompute_dynamics() {
	// Precomputes structure coefficients for 2-D rectangle basis functions.
	// This is performed symbolically, see paper for details

	this.Ck = new SparseMatrix[N];

	// Allocate sparse matrices
	for (int i =0;i<N;i++) {
	    this.Ck[i] = new SparseMatrix(N,N);
	}

	// Calculate the eigenvalues of each basis field.  
	this.eigs = new double[N];
	this.eigs_inv = new double[N];
	this.eigs_inv_root = new double[N];

	for (int i=0;i<N;i++) {
	    int k1 = this.basis_lookup(i,0);
	    int k2 = this.basis_lookup(i,1);
	    this.eigs[i] = (k1*k1 + k2*k2);
	    this.eigs_inv[i] = 1.0/(k1*k1+k2*k2);
	    this.eigs_inv_root[i] = 1.0/Math.sqrt(k1*k1+k2*k2);
	}

	for (int d1 = 0;d1 <N;d1++) {
	    int a1 = this.basis_lookup(d1,0);
	    int a2 = this.basis_lookup(d1,1);
	    int a1_2 = a1*a1;
	    
	    double lambda_a = -(a1*a1 + a2*a2);
	    double inv_lambda_a = -1.0/(a1*a1 + a2*a2);
	    
	    for (int d2 = 0;d2 <N;d2++) {
		
		int b1 = this.basis_lookup(d2,0);
		int b2 = this.basis_lookup(d2,1);

		double lambda_b = -(b1*b1 + b2*b2);
		double inv_lambda_b = -1.0/(b1*b1 + b2*b2);

		int k1 = this.basis_rlookup(a1,a2);
		int k2 = this.basis_rlookup(b1,b2);

		int [][]antipairs = new int[4][2];
		antipairs[0][0] = a1-b1; antipairs[0][1] = a2-b2;
		antipairs[1][0] = a1-b1; antipairs[1][1] = a2+b2;
		antipairs[2][0] = a1+b1; antipairs[2][1] = a2-b2;
		antipairs[3][0] = a1+b1; antipairs[3][1] = a2+b2;
		
		for (int c = 0;c < 4;c++) {
		    int i = antipairs[c][0];
		    int j = antipairs[c][1];

		    int index = this.basis_rlookup(i,j);
		    
		    if (index != -1) {
			double coef = - this.coefdensity(a1,a2,b1,b2,c,0) * inv_lambda_b;
			this.Ck[index].set(k1,k2,coef);
			this.Ck[index].set(k2,k1,coef * -lambda_b/lambda_a);
		    }
		}
	    }
	}
    }
	    
    public double coefdensity(int a1, int b1, int a2,int b2,int c,int tt) {
	if (tt==0) {
            //SS x SS
	    if (c==0) return 0.25 * -(a1*b2 - a2*b1); // --
	    if (c==1) return 0.25 * (a1*b2 + a2*b1);  // -+
	    if (c==2) return 0.25 * -(a1*b2 + a2*b1); // +-
	    if (c==3) return 0.25 * (a1*b2 - a2*b1);  // ++
	} else if (tt==1) {
	    //SC x SS
	    if (c==0) return 0.25 * -(a1*b2 - a2*b1); // --
	    if (c==1) return 0.25 * -(a1*b2 + a2*b1);  // -+
	    if (c==2) return 0.25 * (a1*b2 + a2*b1); // +-
	    if (c==3) return 0.25 * (a1*b2 - a2*b1);  // ++
        } else if (tt==2) {
	    //CS x SS
	    if (c==0) return 0.25 * -(a1*b2 - a2*b1); // --
	    if (c==1) return 0.25 * -(a1*b2 + a2*b1);  // -+
	    if (c==2) return 0.25 * (a1*b2 + a2*b1); // +-
	    if (c==3) return 0.25 * (a1*b2 - a2*b1);  // ++
	} else if (tt==3) {
	    //CS x SS
	    if (c==0) return 0.25 * -(a1*b2 - a2*b1); // --
	    if (c==1) return 0.25 * -(a1*b2 + a2*b1);  // -+
	    if (c==2) return 0.25 * (a1*b2 + a2*b1); // +-
	    if (c==3) return 0.25 * (a1*b2 - a2*b1);  // ++
	} 

	return 0;
    }

    public double[][][] basis_field_2d_rect(int n,int m,double amp) {
	// Calculate Laplacian eigenfunction for eigenvalue (k1,k2) on 2D rectangle

	int a = n;
	int b = m;

	double xfact = 1.0;
	if (n != 0) xfact = -1.0/(a*a + b*b); 
	double yfact = 1.0;
	if (m != 0) yfact = -1.0/(a*a + b*b); 

	double [][][] vf = new double[2][this.mx+1][this.my+1];
	
	double deltax = 3.14f/this.mx;
	double deltay = 3.14f/this.my;
	
	for (int i = 0;i<this.mx+1;i++) {
	    for (int j=0;j<this.my+1;j++) {
		double x = (double)i * deltax;
		double y = (double)j * deltay;

		vf[0][i][j] = amp * xfact * -b * Math.sin(a*x) * Math.cos(b*(y+0.5*deltay));
		vf[1][i][j] = amp * yfact * a * Math.cos(a*(x+0.5*deltax)) * Math.sin(b*y);
	    }

	}

	return vf;
    }

    public double cur_energy() {
	//calculate current energy, sum of squares of coefficients
	// since laplacian eigenfunction basis is orthogonal
	double energy=0.0;
	for (int i=0;i<this.N;i++)
	    energy += this.eigs_inv[i] * (this.coef[i] * this.coef[i]);
	   
	return energy;
    }

    public void set_energy(double desired_e) {
	double cur_e = this.cur_energy();
	double fact = Math.sqrt(desired_e)/Math.sqrt(cur_e);

	for (int i =0;i<this.N;i++) 
	    this.coef[i] *= fact;
    }

    public void attract_particles() {
	// Bunch up particles so we can watch them advect
	for (int i=0;i<this.num_particles;i++) {
	    for (int j=0;j<this.particle_tail_length;j++) {
		particles[i][j][0] *= 0.2;
		particles[i][j][0] += 0.1;
		particles[i][j][1] *= 0.2;
		particles[i][j][1] += 0.4;
	    }
	}
    }

    public void add_particles(int n) {
	this.num_particles = n;
	this.particles=new double[n][this.particle_tail_length][2];
	this.particle_index = new int[n];

	// Seed the particles in random positions
	for (int i=0;i<this.num_particles;i++) {

	    particles[i][0][0] = Math.random();
	    for (int j=0;j<this.particle_tail_length;j++) 
		particles[i][j][0] = particles[i][0][0];

	    particles[i][0][1] = Math.random();
	    for (int j=0;j<this.particle_tail_length;j++) 
		particles[i][j][1] = particles[i][0][1];
	}
    }

    public double getInterpolatedValue(double x, double y, int index) {
	int i = (int)Math.floor(x);
	int j = (int)Math.floor(y);

	double tot = 0.0;
	int den = 0;
	
	if (i >= 0 && i <= mx && j >= 0 && j <= my) {
	    tot += (i+1-x) * (j+1-y) * this.vfield[index][i][j];
	    den++;
	}
	if (i+1 >= 0 && i+1 <= mx && j >= 0 && j <= my) {
	    tot += (x - i) * (j+1-y) * this.vfield[index][i+1][j];
	    den++;
	}

	if (i >= 0 && i <= mx && j+1 >= 0 && j+1 <= my) {
	    tot += (i+1-x) * (y-j) * this.vfield[index][i][j+1];
	    den++;
	}

	if (i+1 >= 0 && i+1 <= mx && j+1 >= 0 && j+1 <= my) {
	    tot += (x-i) * (y-j) * this.vfield[index][i+1][j+1];
	    den++;
	}

	if (den==0) return 0;

	tot = tot/(double)den;

	return tot;
    }

    public double[] vel_at_bilinear(double xx, double yy) {
	double []v = new double[2];
	xx *= this.mx;
	yy *= this.my;
	
	v[0] = getInterpolatedValue(xx,yy-0.5,0);
	v[1] = getInterpolatedValue(xx-0.5,yy,1);
	
	return v;
    }

    public double[] vel_at_cubic(double xx,double yy) {
	double []v = new double[2];
	
	double []d = new double[4];
	double []f = new double[4];
	double tk;
	
	xx *= this.mx;
	yy *= this.my;

	//calculate velocity at x,y
	int k=1;

	int []x = new int[4];
	x[k] = clampi((int)Math.floor(xx),0,mx);
	x[k+1] = clampi(x[k] + 1, 0,mx);
	x[k+2] = clampi(x[k] + 2, 0, mx);
	x[k-1] = clampi(x[k] -1,0,mx);

	int []y = new int[4];
	y[k] = clampi((int)Math.floor(yy),0,my);
	y[k+1] = clampi(y[k] + 1, 0,my);
	y[k+2] = clampi(y[k] + 2, 0, my);
	y[k-1] = clampi(y[k] -1,0,my);

	//x component
	f[k-1] = this.vfield[0][x[k-1]][y[k]];
	f[k] = this.vfield[0][x[k]][y[k]];
	f[k+1] = this.vfield[0][x[k+1]][y[k]];
	f[k+2] = this.vfield[0][x[k+2]][y[k]];
	
	tk = xx - x[k];
	
	v[0] = f[k-1]*(-0.5*tk + tk*tk - 0.5 * tk*tk*tk)
	    +   f[k]*(1.0 - (5.0/2.0)*tk*tk + (3.0/2.0)*tk*tk*tk)
	    +   f[k+1]*(0.5 *tk + 2*tk*tk - (3.0/2.0)*tk*tk*tk)
	    +   f[k+2] *(-0.5 * tk*tk + 0.5*tk*tk*tk);

	//y component
	f[k-1] = this.vfield[1][x[k]][y[k-1]];
	f[k] = this.vfield[1][x[k]][y[k]];
	f[k+1] = this.vfield[1][x[k]][y[k+1]];
	f[k+2] = this.vfield[1][x[k]][y[k+2]];
	
	tk = yy - y[k];
	v[1] = f[k-1]*(-0.5*tk + tk*tk - 0.5 * tk*tk*tk)
	    +   f[k]*(1.0 - (5.0/2.0)*tk*tk + (3.0/2.0)*tk*tk*tk)
	    +   f[k+1]*(0.5 *tk + 2*tk*tk - (3.0/2.0)*tk*tk*tk)
	    +   f[k+2] *(-0.5 * tk*tk + 0.5*tk*tk*tk);

	return v;
    }

    public int clampi(int f, int a, int b) {
	if (f < a) return a;
	if (f > b) return b;
	return f;
    }

    public double clampd(double f, double a, double b) {
	if (f < a) return a;
	if (f > b) return b;
	return f;
    }

    public void advect_density() {
	double [][]density_new = new double[this.dmx][this.dmy];

	double pdt = this.dt * this.pdt_mult;
	boolean RK2 = false;

	for (int i=0;i<this.dmx;i++) {
	    for (int j=0;j<this.dmy;j++) {
		double x = ((double)i +0.5)/this.dmx;
		double y = ((double)j +0.5)/this.dmy;
		
		double nx=0;
		double ny=0;

		if (RK2) {
		    double [] v0 = vel_at_bilinear(x,y);
		    double [] v1 = vel_at_bilinear(x-0.666*pdt*v0[0],y-0.666*pdt*v0[1]);
		    
		    nx = x - pdt*(v0[0] + 3*v1[0])/4.0;
		    ny = y - pdt*(v0[1] + 3*v1[1])/4.0;
		} else {
		    double []v = vel_at_bilinear(x,y);
		    
		    nx = x - pdt * v[0];
		    ny = y - pdt * v[1];
		}

		density_new[i][j] = density_at(nx,ny);
		
	    }
	}

	this.density = density_new;

    }

    public double density_at(double xxx, double yyy) {
	//bilinear interpolation of density
	
	double x  = xxx *this.dmx;
	double y = yyy * this.dmy;

	double xx = clampd(x-0.5,0,(double)(dmx-1));
	double yy = clampd(y-0.5,0,(double)(dmy-1));
	
	int x1 = clampi((int)xx,0,dmx-1);
	int x2 = clampi((int)xx + 1,0,dmx-1);

	int y1 = clampi((int)yy,0,dmy-1);
	int y2 = clampi((int)yy + 1,0,dmy-1);

	double b1 = this.density[x1][y1];
	double b2 = this.density[x2][y1] - this.density[x1][y1];
	double b3 = this.density[x1][y2] - this.density[x1][y1];
	double b4 = this.density[x1][y1] - this.density[x2][y1] - this.density[x1][y2] + this.density[x2][y2];
	
	double dx = xx - (double)x1;
	double dy = yy - (double)y1;

	double tot =  b1 + b2*dx + b3*dy + b4*dx*dy;

	return tot;
    }

    public void advect_particles() {
	//Advects particles using RK4 and cubic velocity interpolation

	double pdt = this.dt * this.pdt_mult;
	
	boolean RK4 = true;
	boolean RK2 = false;
	boolean Euler = false;
	for (int i=0;i<this.num_particles;i++) {
	    int pindex = this.particle_index[i];
	    double x = this.particles[i][pindex][0];
	    double y = this.particles[i][pindex][1];

	    double nx=0;
	    double ny=0;
	    if (RK4) {

		double []v0 = vel_at_bilinear(x,y);
		double []v1 = vel_at_bilinear(x+0.5*pdt*v0[0],y+0.5*pdt*v0[1]);
		double []v2 = vel_at_bilinear(x+0.5*pdt*v1[0],y+0.5*pdt*v1[1]);
		double []v3 = vel_at_bilinear(x+pdt*v2[0],y+pdt*v2[1]);
	    
		nx = x + pdt * (v0[0] + 2*v1[0] +2*v2[0] + v3[0])/6.0;
		ny = y + pdt * (v0[1] + 2*v1[1] +2*v2[1] + v3[1])/6.0;
	    } else if (RK2) {
		double []v0 = vel_at_bilinear(x,y);
		double []v1 = vel_at_bilinear(x-0.666*pdt*v0[0], y-0.666*pdt*v0[1]);

		nx = x + pdt* (v0[0] + 3*v1[0])/4.0;
		ny = y + pdt* (v0[1] + 3*v1[1])/4.0;
	    } else if (Euler) {
		double []v0 = vel_at_bilinear(x,y);
		nx = x + pdt* v0[0];
		ny = y + pdt* v0[1];
	    }

	    nx = clampd(nx,margin,1.0-margin);
	    ny = clampd(ny,margin,1.0-margin);
	    int nextindex = (pindex+1) % this.particle_tail_length;
	    this.particles[i][nextindex][0] = nx;
	    this.particles[i][nextindex][1] = ny;
	    this.particle_index[i] = nextindex;
	}
    }

    public static void main (String[] Argv) {
	System.out.println("Fluid Simulation using Laplacian Eigenfunctions");

	LE fs = new LE();
	
	fs.coef[0]=1.0;
	fs.expand_basis();
	//fs.print_field();
    }

    public void print_field() {
	for (int i=0;i<this.mx;i++) {
	    for (int j=0;j<this.my;j++) {
		System.out.printf("[%d,%d] -> %f\n",i,j, this.vfield[0][i][j]);
	    }
	}
    }

    public int basis_lookup(int index,int component) {
	return this.basis_lookup_table[index][component];
    }

    public int basis_rlookup(int k1,int k2) {
	if (k1 > this.N_sqrt || k1 < 1 || k2 > this.N_sqrt || k2 < 1) {
	    //these fields do not exist
	    return -1;
	}
	return this.basis_rlookup_table[k1][k2];
    }

    public void expand_basis() {
	// Calculate superposition of basis fields
	
	this.vfield = new double[2][mx+1][my+1];

	for (int k=0;k<this.N;k++) {
	    for (int i=0;i<this.mx+1;i++) {
		for (int j=0;j<this.my+1;j++) {
		    this.vfield[0][i][j]  += this.coef[k] * this.vel_basis[k][0][i][j];
		    this.vfield[1][i][j]  += this.coef[k] * this.vel_basis[k][1][i][j];
		}
	    }
	}
    }
    public void fill_lookup_table() {
	// Assume N is a perfect square, and use all basis fields
	// with eigenvalues (k1,k2) up to (N^0.5, N^0.5)

	this.N_sqrt = (int) Math.floor(Math.sqrt((double)N));

	this.basis_lookup_table = new int[this.N][2];
	this.basis_rlookup_table = new int[N_sqrt+1][N_sqrt+1];

	//Initialize lookup table to -1, meaning this (k1,k2) basis field does not exist;
	for (int k1 = 0; k1 <= this.N_sqrt; k1++) {
	    for (int k2 = 0; k2 <= this.N_sqrt; k2++) {
		this.basis_rlookup_table[k1][k2] = -1;
	    }
	}

	int index = 0;
	for (int k1 = 0; k1 <= this.N_sqrt; k1++) {
	    for (int k2 = 0; k2 <= this.N_sqrt; k2++) {
		if (k1 > this.N_sqrt || k1 < 1 || k2 > this.N_sqrt || k2 < 1) {
		    //these fields do not exist
		    continue;
		}

		this.basis_lookup_table[index][0] = k1;
		this.basis_lookup_table[index][1] = k2;

		this.basis_rlookup_table[k1][k2] = index;
		index += 1;
	    }
	}
    }

    public double dot(double []a, double []b) {
	double res=0;
	for (int i=0;i<a.length;i++)
	    res += a[i]*b[i];

	return res;
    }

    public void dump(String s, double []a) {
	System.out.printf("%s [",s);
	for (int i=0;i<a.length;i++)
	    System.out.printf("%f ",a[i]);

	System.out.printf("]\n");
    }

    public void dump(String s, double a) {
	System.out.printf("%s %f\n",s,a);
    }

    public double [] project_forces(double [][]force_path) {
	double []dw = new double[this.N];

	for (int i=0;i<this.N;i++) {
	    double tot = 0;

	    int a = this.basis_lookup(i,0);
	    int b = this.basis_lookup(i,1);

	    double xfact = 1.0;
	    if (a !=0)
		xfact = -1.0/(a*a+b*b);
	    double yfact = 1.0;
	    if (b !=0)
		yfact = -1.0/(a*a+b*b);

	    for (int j=0;j<force_path.length-1;j++) {
		double x = force_path[j][0];
		double y = force_path[j][1];
		double fx = force_path[j][2];
		double fy = force_path[j][3];

		if (x >= 1.00001 || x <= -0.00001 || y >= 1.00001 || y <= -0.00001)
		    continue;

		x *= 3.14;
		y *= 3.14;
	    
		double vx = this.dt * xfact * -b * Math.sin(a*x) * Math.cos(b*y);
		double vy = this.dt * yfact * a * Math.cos(a*x) * Math.sin(b*y);

		tot += (vx*fx + vy*fy);
	    }
	    dw[i] = tot;
	}
	return dw;
    }

    public void stir(double [][]force_path) {
	// Calculate the projected forces, and incorporate them on the next timestep
	double []dw = this.project_forces(force_path);
	for (int i=0;i<this.N;i++)
	    this.forces_dw[i] += dw[i];
    }
}
