#! /usr/bin/python2.7

from time import sleep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#planet1 parameters
a2 = 3**(2/3)  #0.63
e2 = 0.0
w2 = 0.0*np.pi	# Argument of pericentre
f2 = 0.0*np.pi	# Initial true anomaly

#planet2 parameters
a1 = 2**(2/3)  #1.0
e1 = 0.3
w1 = 0.0*np.pi
f1 = 0.0*np.pi	# Initial true anomaly

# Display parameters
tTail = 0.2		# Fraction of a time unit to draw the tail for
resTail = 10			# Display resolution of the tail
tSecond = 0.5		# Time units per second
fps = 50			# Frames per second
tEndSeconds = 60.0	# Wall-clock time to animate up to before reset
resOrbit = 360		# Display resolution of the orbit

# Radius given a, e and true anomaly
def r(a, e, f):
	return (a*(1-e**2)/(1+e*np.cos(f)))

# Return the true anomaly given the mean anomaly
def fFromM(M, e):
	E = EfromM(M, e)
	f = fFromE(E, e)
	return f

def EfromM(M, e, iterations=20):
	# Iteratively solve Kepler's equation
	# M = E - e*sin(E)
	E = M	# E0
	for i in range(1, iterations+1):
		E = M + e*np.sin(E)
	return E

def MfromE(E, e):
	M = E - e*np.sin(E)
	return M

def fFromE(E, e):
	# Solve cos(f) = (cos(E)-e)/(1 - e*cos(E))
	#f = np.arccos((np.cos(E)-e)/(1-e*np.cos(E)))
	#f = np.where(E % (2.0*np.pi) > np.pi,
	#             2.0*np.pi - f,
	#             f)
	f = 2.0*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2.0))
	f = f % (2.0*np.pi)
	return f

def Efromf(f, e):
	E = 2.0*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(f/2.0))
	E = E % (2.0*np.pi)
	return E

# Compute the time at which the last pericentre passage occurred
def pericenterPassageTime(a, e, f0, w0):
	E0 = Efromf(f0, e)
	M0 = MfromE(E0, e)
	# Compute the mean motion
	n = a**(-1.5)
	t0 = -M0 / n
	return t0

# Return the true anomaly given the semi-major axis, initial true anomaly, and time
def fFromt(a, e, t0, t):
	# Computer the mean motion
	n = a**(-1.5)
	# Compute the mean anomaly
	M = n*(t - t0)
	E = EfromM(M, e)
	f = fFromE(E, e)
	return f

# Animation of the orbit
fig, ax = plt.subplots(figsize=(6,6))
ax.axis([-3,3,-3,3])
ax.set_aspect(1.0)


x = np.arange(0, tTail, tTail/resTail)
fullArray = np.arange(0, 2.0*np.pi, 2.0*np.pi/resOrbit)
orbit1, = ax.plot(np.cos(x), np.sin(x), 'k:')
orbit2, = ax.plot(np.cos(x), np.sin(x), 'k:')
planet1, = ax.plot(0.0, 0.0, 'bo')
planet2, = ax.plot(0.0, 0.0, 'ro')
line1, = ax.plot(np.cos(x), np.sin(x), 'b-', linewidth=2)
line2, = ax.plot(np.cos(x), np.sin(x), 'r-', linewidth=2)
star, = ax.plot(0.0, 0.0, 'k*')

def animate(i):
	t = i * tSecond / fps
	orbit1.set_ydata(r(a1, e1, fullArray)*np.sin(fullArray + w1))
	orbit1.set_xdata(r(a1, e1, fullArray)*np.cos(fullArray + w1))
	orbit2.set_ydata(r(a2, e2, fullArray)*np.sin(fullArray + w2))
	orbit2.set_xdata(r(a2, e2, fullArray)*np.cos(fullArray + w2))
	planet1.set_ydata(r(a1, e1, fFromt(a1, e1, t1, t))*np.sin(fFromt(a1, e1, t1, t) + w1))
	planet1.set_xdata(r(a1, e1, fFromt(a1, e1, t1, t))*np.cos(fFromt(a1, e1, t1, t) + w1))
	planet2.set_ydata(r(a2, e2, fFromt(a2, e2, t2, t))*np.sin(fFromt(a2, e2, t2, t) + w2))
	planet2.set_xdata(r(a2, e2, fFromt(a2, e2, t2, t))*np.cos(fFromt(a2, e2, t2, t) + w2))
	line1.set_ydata(r(a1, e1, fFromt(a1, e1, t1, t-x))*np.sin(fFromt(a1, e1, t1, t-x) + w1))
	line1.set_xdata(r(a1, e1, fFromt(a1, e1, t1, t-x))*np.cos(fFromt(a1, e1, t1, t-x) + w1))
	line2.set_ydata(r(a2, e2, fFromt(a2, e2, t2, t-x))*np.sin(fFromt(a2, e2, t2, t-x) + w2))
	line2.set_xdata(r(a2, e2, fFromt(a2, e2, t2, t-x))*np.cos(fFromt(a2, e2, t2, t-x) + w2))
	star.set_ydata(0.0)
	star.set_xdata(0.0)
	return orbit1, orbit2, planet1, planet2, line1, line2, star,

def init():
	orbit1.set_ydata(np.ma.array(x,mask=True))
	orbit1.set_xdata(np.ma.array(x,mask=True))
	orbit2.set_ydata(np.ma.array(x,mask=True))
	orbit2.set_xdata(np.ma.array(x,mask=True))
	planet1.set_ydata(np.ma.array(x,mask=True))
	planet1.set_xdata(np.ma.array(x,mask=True))
	planet2.set_ydata(np.ma.array(x,mask=True))
	planet2.set_xdata(np.ma.array(x,mask=True))
	line1.set_ydata(np.ma.array(x,mask=True))
	line1.set_xdata(np.ma.array(x,mask=True))
	line2.set_ydata(np.ma.array(x,mask=True))
	line2.set_xdata(np.ma.array(x,mask=True))
	star.set_ydata(np.ma.array(x,mask=True))
	star.set_xdata(np.ma.array(x,mask=True))
	return orbit1, orbit2, planet1, planet2, line1, line2, star,

# Compute the times of pericentre passage
t1 = pericenterPassageTime(a1, e1, f1, w1)
t2 = pericenterPassageTime(a2, e2, f2, w2)

ani = animation.FuncAnimation(fig, animate, np.arange(1,fps*tEndSeconds),
       init_func=init, interval=(100.0/fps), blit=True)

plt.show()
FFwriter = animation.FFMpegWriter(fps=50)
ani.save('animation.mp4', writer = FFwriter)