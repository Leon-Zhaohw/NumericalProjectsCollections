#! python3

# Tim Langlois 2016

from tkinter import *
import simpleaudio as sa
from wavefile import WaveWriter, Format
import math
import numpy as np
from scipy.integrate import odeint,ode

# physical constants
CF = 1497.
MU = 8.9e-4
RHO_WATER = 998.
GTH = 1.6e6
GAMMA = 1.4
G = 9.8
SIGMA = 0.072
ETA = 0.84
PATM = 101325.


# a subclass of Canvas for dealing with resizing of windows
class ResizingCanvas(Canvas):
    def __init__(self,parent,**kwargs):
        Canvas.__init__(self,parent,**kwargs)
        self.bind("<Configure>", self.on_resize)
        self.height = self.winfo_reqheight()
        self.width = self.winfo_reqwidth()

    def on_resize(self,event):
        # determine the ratio of old width/height to new width/height
        wscale = float(event.width)/self.width
        hscale = float(event.height)/self.height
        self.width = event.width
        self.height = event.height
        # resize the canvas
        self.config(width=self.width,
                    height=self.height)
        # rescale all the objects tagged with the "all" tag
        self.scale("all",
                   0,
                   0,
                   wscale,
                   hscale)


# A subclass of ResizingCanvas to draw the bubble correctly
class BubbleCanvas(ResizingCanvas):
    def __init__(self,parent,**kwargs):
        ResizingCanvas.__init__(self, parent, **kwargs)
        self.bind("<Configure>", self.on_resize)

        w = self.winfo_reqwidth()
        h = self.winfo_reqheight()

        self.bar = self.create_rectangle(0,
                                         h * .05,
                                         w,
                                         h * .1,
                                         fill="red")

        self.bubble_depth = 2 # in units of bubble diameter
        self.bubble_height = .9 / 50 # in units of window height

        self.draw_bubble()

    def change_depth(self, d):
        self.bubble_depth = float(d)

        self.resize_bubble()

    def change_radius(self, r):
        pass
        #self.bubble_height = 2.0 * float(r)

    def update_height(self, h):
        self.bubble_height = .9 / h
        self.resize_bubble()

    def draw_bubble(self):
        w = self.winfo_reqwidth()
        h = self.winfo_reqheight()

        bub_r = h * self.bubble_height / 2.0
        bub_top = h * .1 + 2.0 * bub_r * self.bubble_depth - bub_r

        self.bubble = self.create_oval(w/2.0 - bub_r,
                                       bub_top,
                                       w/2.0 + bub_r,
                                       bub_top + 2 * bub_r,
                                       outline = 'white',
                                       width = 3)

    def resize_bubble(self):
        w = self.winfo_reqwidth()
        h = self.winfo_reqheight()

        bub_r = h * self.bubble_height / 2.0 # radius
        bub_top = h * .1 + 2.0*bub_r * self.bubble_depth - bub_r

        self.coords(self.bubble,
                    w/2.0 - bub_r,
                    bub_top,
                    w/2.0 + bub_r,
                    bub_top + 2 * bub_r)

    def on_resize(self,event):
        super().on_resize(event)

        # Draw the bubble ourselves to keep it circular
        self.resize_bubble()

    def change_interface(self, t):
        if t == 1:
            self.itemconfig(self.bar, fill="red")
        elif t == 2:
            self.itemconfig(self.bar, fill="blue")


# Class to play bubble sounds
class BubbleSound:
    def __init__(self):
        self.moving_type = 1
        self.interface_type = 1

    def change_type(self, t):
        if t == 1 or t == 2:
            self.moving_type = t

    def change_interface(self, t):
        if t == 1 or t == 2:
            self.interface_type = t

    def bubble_capacitance(self, r, d):
        if self.interface_type == 1: # fluid
            C = r / (1.0 - r / (2 * d) - (r / (2*d))**4)
        else: # Rigid interface
            C = r / (1.0 + r / (2 * d) - (r / (2*d))**4)

        return C


    def minnaert_freq(self, r, d):
        omega = math.sqrt(3 * GAMMA * PATM - 2 * SIGMA * r) / (r * math.sqrt(RHO_WATER))

        f = omega / 2 / math.pi
        return f

    def actual_freq(self, r, d):
        C = self.bubble_capacitance(r, d)

        p0 = PATM
        v0 = 4.0/3.0 * math.pi * r**3
        omega = math.sqrt(4.0 * math.pi * GAMMA * p0 * C / (RHO_WATER * v0))
        f = omega / 2 / math.pi

        return f


    def calc_beta(self, r, w0):
        dr = w0 * r / CF
        dvis = 4. * MU / (RHO_WATER * w0 * r**2)
        phi = 16. * GTH * G / (9 * (GAMMA-1)**2 * w0)
        dth = 2*(math.sqrt(phi - 3.) - (3. * GAMMA - 1.) /
                 (3. * (GAMMA - 1))) / (phi - 4)

        dtotal = dr + dvis + dth

        return w0 * dtotal / math.sqrt(dtotal**2 + 4)

    def jet_forcing(self, r, t):
        cutoff = min(.0006, 0.5 / (3.0 / r))

        if t < 0 or t > cutoff:
            return 0


        jval = (-9 * GAMMA * SIGMA * ETA *
                (PATM + 2 * SIGMA/r) * math.sqrt(1 + ETA**2) /
                (4 * RHO_WATER**2 * r**5 ) * t**2)

        # Convert to radius (instead of fractional radius)
        jval *= r

        # Convert to pressure
        mrp = RHO_WATER * r
        jval *= mrp

        return jval

    # Calculate the bubble terminal velocity according to the paper
    # Rising Velocity for Single Bubbles in Pure Liquids
    # Bax-Rodriguez et al. 2012
    def bubble_terminal_velocity(self, r):
        d = 2 * r

        del_rho = 997. # Density difference between the phases

        # eq 2
        vtpot = 1./36. * del_rho * G * d**2 / MU

        # eq 6
        vt1 = vtpot * math.sqrt(1 + 0.73667 * math.sqrt(G * d) / vtpot)

        # eq 8
        vt2 = math.sqrt(3 * SIGMA / RHO_WATER / d + G * d * del_rho / 2 / RHO_WATER)

        # eq 1
        vt = 1 / math.sqrt( 1 / vt1**2 + 1/vt2**2 )

        return vt

    def bubble_integrator(self, y, t, r, d0, dt, of):
        #[f'; f]

        f = self.jet_forcing(r, t - 0.1)

        d = d0
        if self.moving_type == 2 and t >= 0.1: # rising bubble, calc depth
            vt = self.bubble_terminal_velocity(r)

            d = max(0.51 * 2 * r, d0 - (t - 0.1) * vt)

            #print('vt: ' + str(vt))
            #print('d: ' + str(d))

        # if we let it run too long and the values get very small,
        # the scipy integrator has problems. Might be setting the time step too
        # small? So just exit when the oscillator loses enough energy
        if t > 0.11 and math.sqrt(y[0]**2 + y[1]**2) < 1e-15:
            return [0, 0]

        p0 = PATM + 2.0 * SIGMA / r
        v0 = 4./3. * math.pi * r**3

        w0 = self.actual_freq(r, d) * 2 * math.pi
        k = GAMMA * p0 / v0
        m = k / w0**2

        beta = self.calc_beta(r, w0)

        acc = f / m - 2 * beta * y[0] - w0**2 * y[1]

        if of:
            of.write(str(w0 / 2 / math.pi) + ' ' + str(y[0]) + '\n')

        #print(y)
        if np.isnan(acc) or max(y) > 1e-4:
            print('y: ' + str(y))
            print('f: ' + str(f))
            print('m: ' + str(m))
            print('w0: ' + str(w0))
            print('beta: ' + str(beta))
            print('t: ' + str(t))

            raise Exception('nan')

        return [acc, y[0]]


    def play_bubble(self, r, d, save_file):
        # Integrate the bubble sound into a buffer
        #print('r: ' + str(r))
        #print('d: ' + str(d))
        numsteps = 96001
        dt = 1./(numsteps-1)
        t = np.linspace(0, 1, numsteps-1)
        y0 = [0, 0]

        #of = open('d.txt', 'w')

        sol = odeint(self.bubble_integrator,
                     y0,
                     t,
                     args=(r,d,dt,None),
                     hmax = dt)

        #of.close()
        #print(sol.shape)

        # And play it
        p = sol[:,1]

        #print(p.max())
        #print(p.min())

        p /= max(p.min(), p.max()) * 1.05

        # Save to wave file
        if save_file:
            with WaveWriter('bub.wav', channels=1, samplerate = numsteps-1) as w :
                ps = np.reshape(p, (1, len(p)))
                w.write(ps)

        # Now play the sound
        # Need to convert it first

        p *= 32767
        p = p.astype(np.int16)

        #np.savetxt('p.txt', p)

        # start playback
        play_obj = sa.play_buffer(p, 1, 2, numsteps-1)

        # wait for playback to finish before exiting
        play_obj.wait_done()


# A class to hold all the control buttons/sliders
class ControlPanel:
    def __init__(self, app):
        self.app = app

        # Static and rising buttons
        self.v = IntVar()

        self.static = Radiobutton(app.top,
                                  text="Static",
                                  variable=self.v,
                                  value=1,
                                  command =lambda: app.sound.change_type(1))

        self.static.select()
        self.static.pack(side=LEFT)
        Radiobutton(app.top,
                    text="Rising",
                    variable=self.v,
                    value=2,
                    command =lambda: app.sound.change_type(2)).pack(side=LEFT)

        # Fluid or rigid interface
        self.interface = IntVar()
        self.fluid = Radiobutton(app.top,
                                 text="Fluid",
                                 variable=self.interface,
                                 value=1,
                                 command =lambda: self.change_interface(1))

        self.fluid.select()
        self.fluid.pack(side=LEFT)
        Radiobutton(app.top,
                    text="Rigid",
                    variable=self.interface,
                    value=2,
                    command =lambda: self.change_interface(2)).pack(side=LEFT)

        self.save_bubble = IntVar()
        self.save_bubble_button = Checkbutton(app.top, text="Save Wav", variable=self.save_bubble)
        self.save_bubble_button.pack(side=LEFT)

        self.rise = Button(app.top,
                           text="Play Bubble",
                           command=lambda: app.sound.play_bubble(float(self.radius.get()) / 1000.,
                                                                 float(self.depth.get()) * float(self.radius.get()) / 1000. * 2,
                                                                 int(self.save_bubble.get())))

        self.rise.pack(side=LEFT)

        self.zoom = Scale(app.top,
                          from_=2,
                          to=50,
                          resolution = 0.1,
                          orient=HORIZONTAL,
                          label='Zoom',
                          command=self.change_zoom)

        self.zoom.set(50)

        self.zoom.pack(side=LEFT)

        self.depth = Scale(app.top,
                           from_=0.51,
                           to=50 - .5,
                           resolution = 0.1,
                           orient=HORIZONTAL,
                           label='Depth',
                           command=self.change_depth)

        self.depth.set(1)

        self.depth.pack(side=LEFT)

        self.radius = Scale(app.top,
                            from_=0.5,
                            to=8,
                            resolution = 0.1,
                            orient=HORIZONTAL,
                            label='Radius',
                            command=self.change_radius)

        self.radius.set(1)

        self.radius.pack(side=LEFT)

        self.minnaert_freq = Label(app.top, text = self.minnaert_freq_str(), font = ('Default', 24))
        self.minnaert_freq.pack(side=LEFT)

        self.actual_freq = Label(app.top, text = self.actual_freq_str(), font = ('Helvetica', 24))
        self.actual_freq.pack(side=LEFT)


    def minnaert_freq_str(self):
        #f = 3.0 / (float(self.radius.get()) / 1000.0)
        #return "Minnaert freq: %04.1f Hz  " % f

        r = float(self.radius.get()) / 1000.0
        d = float(self.depth.get()) * r * 2
        f = self.app.sound.minnaert_freq(r, d)

        return "Minnaert freq: %04.1f Hz  " % f

    def actual_freq_str(self):
        r = float(self.radius.get()) / 1000.0
        d = float(self.depth.get()) * r * 2

        return "Actual freq: %04.1f Hz  " % self.app.sound.actual_freq(r, d)

    def change_zoom(self, z):
        h = float(z)
        self.depth.config(to=h - .5)
        self.change_depth(float(self.depth.get()))
        self.app.bubble_canvas.update_height(h)

    def change_interface(self, t):
        self.app.sound.change_interface(t)
        self.app.bubble_canvas.change_interface(t)
        self.minnaert_freq.config(text=self.minnaert_freq_str())
        self.actual_freq.config(text=self.actual_freq_str())

    def change_depth(self, d):
        self.app.bubble_canvas.change_depth(float(d))
        self.minnaert_freq.config(text=self.minnaert_freq_str())
        self.actual_freq.config(text=self.actual_freq_str())

    def change_radius(self, r):
        self.app.bubble_canvas.change_radius(float(r))
        self.minnaert_freq.config(text=self.minnaert_freq_str())
        self.actual_freq.config(text=self.actual_freq_str())


# The main application
class App:
    def __init__(self, master):
        self.m = PanedWindow(master, orient=VERTICAL)
        self.m.pack(fill=BOTH, expand=YES)

        self.top = Frame(self.m)
        self.m.add(self.top)

        self.bubble_canvas = BubbleCanvas(self.m,
                                          width=400,
                                          height=800,
                                          bg='black')

        self.bubble_canvas.pack(fill=BOTH, expand=YES)
        self.m.add(self.bubble_canvas)

        self.sound = BubbleSound()

        self.control_panel = ControlPanel(self)


root = Tk()

app = App(root)

root.mainloop()


