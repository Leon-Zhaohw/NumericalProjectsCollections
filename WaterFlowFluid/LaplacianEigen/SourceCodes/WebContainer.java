/* Tyler de Witt
2012-08-03
*/

import java.applet.*;
import java.awt.event.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import javax.swing.Timer;
import java.util.*;

public class WebContainer extends Applet
    implements MouseListener, MouseMotionListener, KeyListener,  AdjustmentListener, ActionListener
{

    // frame dimensions (dxd pixels)
    int d = 400;

    // The fluid solver 
    LE fs;

    // Animation timer event 
    Timer clock;

    // Controls
    Scrollbar []control_slider;
    TextField []control_text;
    String []control_label;
    Button []control_button;

    // mouse position
    int x, xOld;
    int y, yOld;

    // cell index
    int i, j;

    // cell dimensions
    int dg, dg_2;

    // cell position
    int dx, dy;

    int mouse_button;

    BufferedImage bi;
    Graphics2D big;

    //External forces

    // Mode 1 saves up the forces from a complete mouse drag and projects them in one shot
    // mode 2 projects the forces as the mouse is dragged
    int force_mode = 2;

    // Strength factor of mouse forces;
    double force_mult = 250.0;

    //force vector on mouse drag
    double []p1 = new double[2];
    double []p2 = new double[2];

    ArrayList<Integer> drag_path_x;
    ArrayList<Integer> drag_path_y;

    // Default parameters
    double dt = 0.1;
    int mx = 64;
    int dmx = 64;
    int num_particles=1000;
    int N = 64;
    int N_sqrt = 8;
    int visc_coarse = 0;
    int visc_fine = 100;
    double visc = 0.0;

    boolean show_vfield = false;
    boolean show_particles = true;
    boolean show_density  = true;

    public void reset() {
	this.fs =new LE(this.mx,this.N,this.dmx);

	this.fs.visc = this.visc;
	this.fs.dt = this.dt;
	fs.particle_tail_length=3;
	fs.add_particles(this.num_particles);	
	fs.expand_basis();
    }
    
    public void init()
    {
	this.visc = ((double)this.visc_coarse + (double)this.visc_fine/200.0)/500.0;

	// Reset fluid simulation
	reset();
	
        addMouseMotionListener(this);
        addMouseListener(this);
        addKeyListener(this);
        setFocusable(true);
        bi = (BufferedImage) createImage(d, d);
        big = bi.createGraphics();
	
	this.setLayout(null);

	//UI
	Panel f = new Panel();

	int rows=10;
	int cols=3;

	f.setLayout(null);
	f.setBounds(400,0,300,rows*30);
	
	control_slider = new Scrollbar[10];
	control_text = new TextField[10];
	control_label = new String[10];
	control_button = new Button[10];

	int ci=0;

	// Viscosity Control
	control_label[0]="Viscosity Coarse";
	control_slider[0] = new Scrollbar(Scrollbar.HORIZONTAL,this.visc_coarse,50,0,500);
	control_slider[0].addAdjustmentListener(this);
	control_slider[0].setBounds(100,ci,150,20);


	control_text[0] = new TextField(4);
	control_text[0].setText(String.valueOf(this.visc*1000));
	control_text[0].setSize(40,40);
	control_text[0].setBounds(250,ci,50,20);


	Label l = new Label(control_label[0]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[0]);
	f.add(control_text[0]);


	// Viscosity Control
	ci += 22;
	control_label[4]="Viscosity Fine";
	control_slider[4] = new Scrollbar(Scrollbar.HORIZONTAL,this.visc_fine,50,0,1000);
	control_slider[4].addAdjustmentListener(this);
	control_slider[4].setBounds(100,ci,150,20);


	l = new Label(control_label[4]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[4]);



	// Num Particles Control
	ci += 22;
	control_label[1]="Num Particles";
	control_slider[1] = new Scrollbar(Scrollbar.HORIZONTAL,this.num_particles,1,0,10000);
	control_slider[1].addAdjustmentListener(this);
	control_slider[1].setBounds(100,ci,150,20);

	control_text[1] = new TextField(4);
	control_text[1].setText(String.valueOf(this.num_particles));
	control_text[1].setSize(40,40);
	control_text[1].setBounds(250,ci,50,20);

	l = new Label(control_label[1]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[1]);
	f.add(control_text[1]);


	// Num Basis Fields
	ci += 22;
	control_label[2]="Basis Dimension: N";
	control_slider[2] = new Scrollbar(Scrollbar.HORIZONTAL,this.N_sqrt,1,0,16);
	control_slider[2].addAdjustmentListener(this);
	control_slider[2].setBounds(100,ci,150,20);

	control_text[2] = new TextField(4);
	control_text[2].setText(String.valueOf(this.N));
	control_text[2].setSize(40,40);
	control_text[2].setBounds(250,ci,50,20);

	l = new Label(control_label[2]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[2]);
	f.add(control_text[2]);

	// Grid Resolution
	ci += 22;
	control_label[3]="Velocity Grid Res";
	control_slider[3] = new Scrollbar(Scrollbar.HORIZONTAL,this.mx,1,0,64);
	control_slider[3].addAdjustmentListener(this);
	control_slider[3].setBounds(100,ci,150,20);

	control_text[3] = new TextField(4);
	control_text[3].setText(String.valueOf(this.mx));
	control_text[3].setSize(40,40);
	control_text[3].setBounds(250,ci,50,20);

	l = new Label(control_label[3]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[3]);
	f.add(control_text[3]);


	// Density Grid Resolution
	ci += 22;
	control_label[5]="Density Grid Res";
	control_slider[5] = new Scrollbar(Scrollbar.HORIZONTAL,this.dmx,1,0,200);
	control_slider[5].addAdjustmentListener(this);
	control_slider[5].setBounds(100,ci,150,20);

	control_text[5] = new TextField(4);
	control_text[5].setText(String.valueOf(this.dmx));
	control_text[5].setSize(40,40);
	control_text[5].setBounds(250,ci,50,20);

	l = new Label(control_label[5]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[5]);
	f.add(control_text[5]);




	// Time step size
	ci += 22;
	control_label[6]="Time step";
	control_slider[6] = new Scrollbar(Scrollbar.HORIZONTAL,(int)(this.dt*4000),1,0,1000);
	control_slider[6].addAdjustmentListener(this);
	control_slider[6].setBounds(100,ci,150,20);

	control_text[6] = new TextField(4);
	control_text[6].setText(String.valueOf(this.dt));
	control_text[6].setSize(40,40);
	control_text[6].setBounds(250,ci,50,20);

	l = new Label(control_label[6]);
	l.setBounds(0,ci,100,20);
	f.add(l);
	f.add(control_slider[6]);
	f.add(control_text[6]);

	// Reset Button
	ci += 22;
	control_button[2] = new Button("Reset");
	control_button[2].addActionListener(this);
	control_button[2].setBounds(0,ci,250,20);
	f.add(control_button[2]);


	// Reseed Particles Locally Button
	ci += 22;
	control_button[0] = new Button("Bunch up Particles");
	control_button[0].addActionListener(this);
	control_button[0].setBounds(0,ci,250,20);
	f.add(control_button[0]);

	// Reseed Particles Globally Button
	ci += 22;
	control_button[1] = new Button("Reseed Particles Globally");
	control_button[1].addActionListener(this);
	control_button[1].setBounds(0,ci,250,20);
	f.add(control_button[1]);

	// Show Particles Button
	ci += 22;
	control_button[3] = new Button("Toggle Show Particles");
	control_button[3].addActionListener(this);
	control_button[3].setBounds(0,ci,250,20);
	f.add(control_button[3]);

	// Show Velocity Field Button
	ci += 22;
	control_button[4] = new Button("Toggle Show Velocity Field");
	control_button[4].addActionListener(this);
	control_button[4].setBounds(0,ci,250,20);
	f.add(control_button[4]);

	// Show Velocity Field Button
	ci += 22;
	control_button[5] = new Button("Toggle Show Density Field");
	control_button[5].addActionListener(this);
	control_button[5].setBounds(0,ci,250,20);
	f.add(control_button[5]);

	this.add(f);
    }

    public void actionPerformed(ActionEvent e) {

	if (e.getSource()==this.clock) {
	    long el = System.nanoTime();
	    fs.step();
	    long elapsed_step = (System.nanoTime() - el)/1000;
	    
	    el = System.nanoTime();
	    if (this.show_particles)
		fs.advect_particles();
	    if (this.show_density)
		fs.advect_density();
	    long elapsed_particles = (System.nanoTime() - el)/1000;
	    //System.out.printf("Elapsed microseconds: %d, Step: %d, Particles %d\n",elapsed_step + elapsed_particles, elapsed_step, elapsed_particles);

            repaint();
	} else if (e.getSource()==this.control_button[0]) {
	    //Reseed particles locally
	    fs.attract_particles();
	} else if (e.getSource()==this.control_button[1]) {
	    //Reseed particles globally
	    fs.add_particles(this.num_particles);
	} else if (e.getSource()==this.control_button[2]) {
	    //Reset fluid solver
	    this.reset();
	} else if (e.getSource()==this.control_button[3]) {
	    this.show_particles = !this.show_particles;
	} else if (e.getSource()==this.control_button[4]) {
	    this.show_vfield = !this.show_vfield;
	} else if (e.getSource()==this.control_button[5]) {
	    this.show_density = !this.show_density;
	}
	
    }

    public void adjustmentValueChanged(AdjustmentEvent e) {
	if (e.getSource()==this.control_slider[0] || e.getSource()==this.control_slider[4]) {
	    // Viscosity
	    this.visc_coarse = this.control_slider[0].getValue();
	    this.visc_fine = this.control_slider[4].getValue();
	    this.visc = ((double)this.visc_coarse + (double)this.visc_fine/200.0)/500.0;

	    control_text[0].setText(String.valueOf(this.visc*1000));

	    // update viscosity immediately
	    this.fs.visc = this.visc;

	} else 	if (e.getSource()==this.control_slider[1]) {
	    // Number of particles
	    int val = this.control_slider[1].getValue();
	    this.num_particles = val;
	    control_text[1].setText(String.valueOf(this.num_particles));
	} else 	if (e.getSource()==this.control_slider[2]) {
	    // Basis Dimension
	    int val = this.control_slider[2].getValue();
	    this.N = val*val;
	    this.N_sqrt= val;
	    control_text[2].setText(String.valueOf(this.N));
	} else 	if (e.getSource()==this.control_slider[3]) {
	    // Vel Grid Res
	    int val = this.control_slider[3].getValue();
	    this.mx = val;
	    control_text[3].setText(String.valueOf(this.mx));
	} else 	if (e.getSource()==this.control_slider[5]) {
	    // Density Grid Res
	    int val = this.control_slider[5].getValue();
	    this.dmx = val;
	    control_text[5].setText(String.valueOf(this.dmx));
	} else 	if (e.getSource()==this.control_slider[6]) {
	    // Time step
	    int val = this.control_slider[6].getValue();
	    this.dt = (double)val/4000.0;
	    this.fs.dt = this.dt;
	    control_text[6].setText(String.valueOf(this.dt));
	}
    }

    public void start()
    {
        // calculate cell deimensions
        dg   = d  / fs.mx;
        dg_2 = dg / 2;

	this.clock = new Timer(40,this);
	this.clock.start();
    }


    public void update(Graphics g){ paint(g); }

    public void plot_particles_fancy() {
	big.setColor(Color.black);

	int col_delta  = 255/(fs.particle_tail_length-1);


	for (int i=0;i<fs.num_particles;i++) {
	    Color col = new Color(255,255,255);
	    int pindex = fs.particle_index[i]+1;
	    if (pindex >= fs.particle_tail_length) pindex=0;

	    int count=0;
	    while (count < fs.particle_tail_length-1) {
		int pindex2 = pindex + 1;
		if (pindex2 >= fs.particle_tail_length) pindex2=0;
		
		double x =fs.particles[i][pindex][0];
		double y =fs.particles[i][pindex][1];		
		double x2 =fs.particles[i][pindex2][0];
		double y2 =fs.particles[i][pindex2][1];		

		Color col2 = new Color(col.getRed() - col_delta,
				       col.getGreen() - col_delta,
				       col.getBlue() - col_delta);

		big.setPaint(new GradientPaint((int)(x*d),(int)(y*d),col,(int)(x2*d), (int)(y2*d),col2,true));
		big.drawLine((int)(x*d),(int)(y*d),(int)(x2*d), (int)(y2*d));
		pindex = pindex2;
		col = col2;
		count++;
	    }
	}
		
	big.setPaint(Color.black);
    }

    public void plot_density() {
        // clear screen
        big.setColor(Color.white);
        big.fillRect(0, 0, d, d);
	int c;

	double delta = 1.0/this.fs.dmx;

	for (int i=0;i<this.fs.dmx;i++) {
	    for (int j=0;j<this.fs.dmy;j++) {
		double x = (double)i * delta;
		double y = (double)j * delta;

		c = (int)((1.0 - fs.density[i][j]) * 255);
		if (c < 0) c=0;
		big.setColor(new Color(c,c,c));
		big.fillRect((int)(x*d),(int)(y*d),(int)((x+delta)*d),(int)((y+delta)*d));
	    }
	}
    }

    public void plot_particles() {
	big.setColor(Color.black);

	for (int i=0;i<fs.num_particles;i++) {
	    int pindex = fs.particle_index[i]+1;
	    if (pindex >= fs.particle_tail_length) pindex=0;

	    int count=0;
	    while (count < fs.particle_tail_length-1) {
		int pindex2 = pindex + 1;
		if (pindex2 >= fs.particle_tail_length) pindex2=0;
		
		double x =fs.particles[i][pindex][0];
		double y =fs.particles[i][pindex][1];		
		double x2 =fs.particles[i][pindex2][0];
		double y2 =fs.particles[i][pindex2][1];		

		big.drawLine((int)(x*d),(int)(y*d),(int)(x2*d), (int)(y2*d));
		pindex = pindex2;
		count++;
	    }
	}
		
	big.setPaint(Color.black);
    }

    public void plot_vfield() {
	big.setColor(Color.red);

        for (int i = 1; i <= fs.mx; i++)
        {
            // x position of current cell
            dx = (int) ( (double)(i - 0.5)/fs.mx * d);

            for (int j = 1; j <= fs.my; j++)
            {
                // y position of current cell
		dy = (int) ( (double)(j - 0.5)/fs.my * d);

		if (i % 2 == 1 && j % 2 == 1) {
		    int u = (int)( 30 * fs.vfield[0][i][j]);
		    int v = (int)( 30 * fs.vfield[1][i][j]);
		    big.drawLine(dx, dy, dx+u, dy+v);
		}
		
	    }
	}
    }
    
    public void paint(Graphics g)
    {
	Graphics2D g2 = (Graphics2D) g;

	// clear screen
	big.setColor(Color.white);
	big.fillRect(0, 0, d, d);

	if (this.show_density)
	    this.plot_density();
	if (this.show_vfield)
	    this.plot_vfield();
	if (this.show_particles)
	    this.plot_particles();


	g2.drawImage(bi, null, 0, 0);
    }

    public void keyPressed(KeyEvent e)
    {
    }

    public void mousePressed(MouseEvent e)
    {
	mouse_button = e.getButton();

	if (this.force_mode==1) {
	    // For forces mode 1, save up the forces from the complete mouse drag, and project in one shot
	    this.drag_path_x = new ArrayList<Integer>();
	    this.drag_path_y = new ArrayList<Integer>();
	    drag_path_x.add(new Integer(e.getX()));
	    drag_path_y.add(new Integer(e.getY()));
	} else {
	    // For forces mode 2, project forces as mouse is dragged
	    this.p1[0] = (double)e.getX() / (double)this.d;
	    this.p1[1] = (double)e.getY() / (double)this.d;
	}
    }

    public void mouseDragged(MouseEvent e)
    {
	if (mouse_button == 3) {
	    //add density
	    int i = (e.getX()* this.fs.dmx) / this.d;
	    int j = (e.getY()* this.fs.dmy) / this.d;

	    if (i >=0 && i < this.fs.dmx && j >=0 && j < this.fs.dmy) {
		this.fs.density[i][j]=100;
	    }

	}  else {
	    if (this.force_mode==1) {
		// For forces mode 1, save up the forces from the complete mouse drag, and project in one shot
		drag_path_x.add(new Integer(e.getX()));
		drag_path_y.add(new Integer(e.getY()));
	    } else {

		// For forces mode 2, project forces as the mouse is dragged
		this.p2[0] = (double)e.getX() /(double)this.d;
		this.p2[1] = (double)e.getY() /(double)this.d;

		double [][]force_path = new double[2][4];
		force_path[0][0] = p1[0];
		force_path[0][1] = p1[1];
		force_path[0][2] = (p2[0] - p1[0]) * this.force_mult;
		force_path[0][3] = (p2[1] - p1[1]) * this.force_mult;

		fs.stir(force_path);

		this.p1[0] = p2[0];
		this.p1[1] = p2[1];
	    }
	}
    }

    public void updateLocation(MouseEvent e)
    {
    }

    // fulfill mouse interface requirements
    public void mouseReleased(MouseEvent e){
	if (force_mode==1) {
	    // For forces mode 1, save up the forces from the complete mouse drag, and project in one shot

	    double [][]force_path= new double[drag_path_x.size()][4];

	    Iterator<Integer> itr_x=drag_path_x.iterator();
	    Iterator<Integer> itr_y=drag_path_y.iterator();

	    int i =0;
	    while(itr_x.hasNext()) {
		Integer x = itr_x.next();
		Integer y = itr_y.next();
		force_path[i][0] = x.doubleValue() / (double)this.d;
		force_path[i][1] = y.doubleValue() / (double)this.d;
		i++;
	    }

	    for (i=0;i<force_path.length-1;i++) {
		force_path[i][2] = (force_path[i+1][0] - force_path[i][0]) * this.force_mult;
		force_path[i][3] = (force_path[i+1][1] - force_path[i][1]) * this.force_mult;
	    }

	    //apply forces to simulation
	    fs.stir(force_path);
	}
    }

    public void mouseMoved(MouseEvent e){}
    public void mouseClicked(MouseEvent e){}
    public void mouseExited(MouseEvent e){}
    public void mouseEntered(MouseEvent e){}

    // fulfill key interface requirements
    public void keyTyped(KeyEvent e){}
    public void keyReleased(KeyEvent e){}
}
