"""
Part 1)
    This cell will create the appropriate lengths for each pendulum
    based on the length of the first pendulum and store them in an 
    array in order to get the animation to work correctly. We use the
    equation derived in the cell above.
"""
from scipy.constants import g # get the gravitational force
from math import pi, pow, sqrt
from numpy import zeros, arange, sin, cos, array

# Initial values:
T_max = 2*60 # amount of time in seconds until repeating pattern
n = 1; # the pendulum we are at (first one here)
L = 10 # length in meters of the very first pendulum // default: 10
k = T_max/(2*pi)*sqrt(g/L)-n-1 # constant based on the length of the pendulum

# Create the array of pendulum lengths:
pendulum_length = zeros([16],float); # array of pendulum lengths with 16 pendulums
pendulum_length[0] = L # assign the first pendulum with initial length L

# Calculates the correct lengths of each next pendulum:
index = 1 # next pendulum stored in the array
for n in arange(2,17): # there are 15 pendulums to calculate the correct lengths
    pendulum_length[index] = g*pow(T_max/(2*pi*(k+n+1)), 2) # length of the pendulum in meters
    index += 1 # increnment the index to store value of next pendulum length in the array
    
"""
This cell displays the lengths of each pendulum on the apparatus.
"""
print("Pendulum Lengths")
# run through pendulum_length list
pend = 1 # the pendulum we are at
for i in pendulum_length:
    length = f'{i:.2f}' # making each length have 2 significant figures
    print("\tPendulum {}: {} meters" .format(pend,length))
    pend += 1 # go to the next pendulum
    
"""
Part 2)
    This cell is where I created each pendulum individually.
    There are 16 total pendulums that are created here. Each
    pendulum alternates in color. Every other pendulum is white.
"""
from vpython import *

# Set the screen width and height:
scene = canvas(); scene.width = 800; scene.height = 400

# Making the pendulum holder:
Length = 30 # length of holder
holder = cylinder(pos=vector(0,0,0), axis=vector(0,0,-Length), radius=1) # holder position

ball_1 = sphere(pos=vector(0,-pendulum_length[0],0), radius=1, color=color.orange) 
string_1 = cylinder(pos=vector(0,0,0), axis=vector(0,-pendulum_length[0],0), radius=0.06)

ball_2 = sphere(pos=vector(0,-pendulum_length[1],-2), radius=1, color=color.white)
string_2 = cylinder(pos=vector(0,0,-2), axis=vector(0,-pendulum_length[1],0), radius=0.06)

ball_3 = sphere(pos=vector(0,-pendulum_length[2],-4), radius=1, color=color.blue)
string_3 = cylinder(pos=vector(0,0,-4), axis=vector(0,-pendulum_length[2],0), radius=0.06)

ball_4 = sphere(pos=vector(0,-pendulum_length[3],-6), radius=1, color=color.white)
string_4 = cylinder(pos=vector(0,0,-6), axis=vector(0,-pendulum_length[3],0), radius=0.06)

ball_5 = sphere(pos=vector(0,-pendulum_length[4],-8), radius=1, color=color.green)
string_5 = cylinder(pos=vector(0,0,-8), axis=vector(0,-pendulum_length[4],0), radius=0.06)

ball_6 = sphere(pos=vector(0,-pendulum_length[5],-10), radius=1, color=color.white)
string_6 = cylinder(pos=vector(0,0,-10), axis=vector(0,-pendulum_length[5],0), radius=0.06)

ball_7 = sphere(pos=vector(0,-pendulum_length[6],-12), radius=1, color=color.red)
string_7 = cylinder(pos=vector(0,0,-12), axis=vector(0,-pendulum_length[6],0), radius=0.06)

ball_8 = sphere(pos=vector(0,-pendulum_length[7],-14), radius=1, color=color.white)
string_8 = cylinder(pos=vector(0,0,-14), axis=vector(0,-pendulum_length[7],0), radius=0.06)

ball_9 = sphere(pos=vector(0,-pendulum_length[8],-16), radius=1, color=color.purple)
string_9 = cylinder(pos=vector(0,0,-16), axis=vector(0,-pendulum_length[8],0), radius=0.06)

ball_10 = sphere(pos=vector(0,-pendulum_length[9],-18), radius=1, color=color.white)
string_10 = cylinder(pos=vector(0,0,-18), axis=vector(0,-pendulum_length[9],0), radius=0.06)

ball_11 = sphere(pos=vector(0,-pendulum_length[10],-20), radius=1, color=color.orange)
string_11 = cylinder(pos=vector(0,0,-20), axis=vector(0,-pendulum_length[10],0), radius=0.06)

ball_12 = sphere(pos=vector(0,-pendulum_length[11],-22), radius=1, color=color.white)
string_12 = cylinder(pos=vector(0,0,-22), axis=vector(0,-pendulum_length[11],0), radius=0.06)

ball_13 = sphere(pos=vector(0,-pendulum_length[12],-24), radius=1, color=color.blue)
string_13 = cylinder(pos=vector(0,0,-24), axis=vector(0,-pendulum_length[12],0), radius=0.06)

ball_14 = sphere(pos=vector(0,-pendulum_length[13],-26), radius=1, color=color.white)
string_14 = cylinder(pos=vector(0,0,-26), axis=vector(0,-pendulum_length[13],0), radius=0.06)

ball_15 = sphere(pos=vector(0,-pendulum_length[14],-28), radius=1, color=color.green)
string_15 = cylinder(pos=vector(0,0,-28), axis=vector(0,-pendulum_length[14],0), radius=0.06)

ball_16 = sphere(pos=vector(0,-pendulum_length[15],-30), radius=1, color=color.white)
string_16 = cylinder(pos=vector(0,0,-30), axis=vector(0,-pendulum_length[15],0), radius=0.06)

"""
Part 3)
    This cell holds the movement of each individual pendulum 
    in a movement function. Each function contains the length
    of that particular pendulum and also the angle for which 
    that pendulum will sway based on its length. There is a
    damping force inside each function so that the pendulums
    don't swing forever and thus creating a more realistic
    animation.
"""

def movement_of_pendulum_1(t,mu,theta_o):
    L = pendulum_length[0]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t) # equation of moving pendulum
    
    ball_1.pos = vector(L*sin(theta), -L*cos(theta), 0) # Update the ball's position
    string_1.axis = vector(L*sin(theta), -L*cos(theta), 0) # updating the string with the ball
    
def movement_of_pendulum_2(t,mu,theta_o):     
    L = pendulum_length[1]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)

    ball_2.pos = vector(L*sin(theta), -L*cos(theta), -2) 
    string_2.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_3(t,mu,theta_o):   
    L = pendulum_length[2]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
 
    ball_3.pos = vector(L*sin(theta), -L*cos(theta), -4) 
    string_3.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_4(t,mu,theta_o):  
    L = pendulum_length[3]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
    
    ball_4.pos = vector(L*sin(theta), -L*cos(theta), -6) 
    string_4.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_5(t,mu,theta_o):  
    L = pendulum_length[4]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
        
    ball_5.pos = vector(L*sin(theta), -L*cos(theta), -8)
    string_5.axis = vector(L*sin(theta), -L*cos(theta), 0)
    
def movement_of_pendulum_6(t,mu,theta_o):
    L = pendulum_length[5]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
    
    ball_6.pos = vector(L*sin(theta), -L*cos(theta), -10)
    string_6.axis = vector(L*sin(theta), -L*cos(theta), 0)
    
def movement_of_pendulum_7(t,mu,theta_o): 
    L = pendulum_length[6]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)

    ball_7.pos = vector(L*sin(theta), -L*cos(theta), -12)
    string_7.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_8(t,mu,theta_o):    
    L = pendulum_length[7]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
 
    ball_8.pos = vector(L*sin(theta), -L*cos(theta), -14)
    string_8.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_9(t,mu,theta_o):  
    L = pendulum_length[8]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
    
    ball_9.pos = vector(L*sin(theta), -L*cos(theta), -16)
    string_9.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_10(t,mu,theta_o):  
    L = pendulum_length[9]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
        
    ball_10.pos = vector(L*sin(theta), -L*cos(theta), -18)
    string_10.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_11(t,mu,theta_o):
    L = pendulum_length[10]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)

    ball_11.pos = vector(L*sin(theta), -L*cos(theta), -20)
    string_11.axis = vector(L*sin(theta), -L*cos(theta), 0)
    
def movement_of_pendulum_12(t,mu,theta_o):   
    L = pendulum_length[11]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)

    ball_12.pos = vector(L*sin(theta), -L*cos(theta), -22)
    string_12.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_13(t,mu,theta_o):    
    L = pendulum_length[12]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
 
    ball_13.pos = vector(L*sin(theta), -L*cos(theta), -24)
    string_13.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_14(t,mu,theta_o):  
    L = pendulum_length[13]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
    
    ball_14.pos = vector(L*sin(theta), -L*cos(theta), -26)
    string_14.axis = vector(L*sin(theta), -L*cos(theta), 0)

def movement_of_pendulum_15(t,mu,theta_o):  
    L = pendulum_length[14]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
        
    ball_15.pos = vector(L*sin(theta), -L*cos(theta), -28)
    string_15.axis = vector(L*sin(theta), -L*cos(theta), 0)
    
def movement_of_pendulum_16(t,mu,theta_o):  
    L = pendulum_length[15]; theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)
        
    ball_16.pos = vector(L*sin(theta), -L*cos(theta), -30)
    string_16.axis = vector(L*sin(theta), -L*cos(theta), 0)
    
"""
Part 4)
    This is the cell that creates the time interval of the pendulums.
    The for loop calls the pendulum functions to have them all swing
    at the same time together to complete the animation.
"""
# Time interval and step of swing in seconds
t1 = 0; t2 = 2*60; tstep = 0.01
tpoints = arange(t1,t2,tstep)

# The damping constant
mu = 0.02 # default 0.02

# Setting the pudulum constants:
theta_o = 1.2 # angle in radians

"""
The below animation will run below the cell that contains (Part 2) above. 
To see the pendulum wave machine, scrole up to the visualization 
after running this cell.
"""
for t in tpoints: # run through the time interval
    rate(275)
    movement_of_pendulum_1(t,mu,theta_o) 
    movement_of_pendulum_2(t,mu,theta_o) 
    movement_of_pendulum_3(t,mu,theta_o) 
    movement_of_pendulum_4(t,mu,theta_o) 
    movement_of_pendulum_5(t,mu,theta_o) 
    movement_of_pendulum_6(t,mu,theta_o) 
    movement_of_pendulum_7(t,mu,theta_o) 
    movement_of_pendulum_8(t,mu,theta_o) 
    movement_of_pendulum_9(t,mu,theta_o) 
    movement_of_pendulum_10(t,mu,theta_o) 
    movement_of_pendulum_11(t,mu,theta_o) 
    movement_of_pendulum_12(t,mu,theta_o) 
    movement_of_pendulum_13(t,mu,theta_o)
    movement_of_pendulum_14(t,mu,theta_o) 
    movement_of_pendulum_15(t,mu,theta_o)
    movement_of_pendulum_16(t,mu,theta_o)
    
"""
Part 5)
    This cell shows plots of the pendulums using the Euler method, 
    the second-order Runge-Kutta method, and the fourth-order Runge-Kutta 
    method that was in chapter 8 of our computational physics textbook.
"""
import matplotlib.pyplot as plt

# The pendulum constants determined previously:
L = 10 # pendulum length in meters // default is 10
N = (t2-t1)/tstep # number of sample points
h = (t2-t1)/N # the step size
mu = 0.02 # default 0.02
theta_o = 1.2 # angle in radians of release
theta = theta_o*exp(-mu*t)*cos(sqrt(g/L)*t)

def f(r,t): # function to solve the ODE 
    theta = r[0]
    omega = r[1]
    ftheta = omega
    fomega = -(g/L)*sin(theta)
    return array([ftheta, fomega],float)

#############################
###### Euler's Method #######
#############################
theta_points1 = [] # list to hold theta values
r1 = array([theta_o,0],float) # aray values to be stored in theta_points list for plot
for t in tpoints: # Euler's method
    theta_points1.append(r1[0])
    r1 += h*f(r1,t)

# Visualizing the plot with Euler's method:
plt.title("Pendulum Wave using Euler's Method") 
plt.xlabel("Time in Seconds"); plt.ylabel("Theta Points in Radians")
plt.plot(tpoints, theta_points1, 'b')
plt.grid(True); plt.show()

#########################################
#### Second-Order Runge-Kutta Method ####
#########################################
theta_points2 = [] # list to hold theta values
r2 = array([theta_o,0],float) # aray values to be stored in theta_points list for plot
for t in tpoints: # 2nd-order Runge-Kutta method
    theta_points2.append(r2[0])
    k1 = h*f(r2,t)
    k2 = h*f(r2+0.5*k1,t+0.5*h)
    r2 += k2
    
# Visualizing the plot with Second-Order Runge-Kutta method:
plt.plot(tpoints, theta_points2, 'g')
plt.title("Pendulum Wave Using Second-Order Runge-Kutta Method")
plt.xlabel("Time in Seconds"); plt.ylabel("Theta Points in Radians")
plt.grid(True); plt.show()

#########################################
#### Fourth-Order Runge-Kutta Method ####
#########################################
theta_points3 = [] # empty list of theta values
r3 = array([theta_o,0],float) # aray values to be stored in theta_points list for plot
for t in tpoints: # 4th-order Runge-Kutta
    theta_points3.append(r3[0])
    k1 = h * f(r3,t)
    k2 = h * f(r3+0.5*k1,t+0.5*h)
    k3 = h * f(r3+0.5*k2,t+0.5*h)
    k4 = h * f(r3+k3,t+h)
    r3 += (1/6)*(k1+2*k2+2*k3+k4)

# Visualizing the plot with Fourth-Order Runge-Kutta method:
plt.plot(tpoints, theta_points3, 'r')
plt.title("Pendulum Wave Using Fourth-Order Runge-Kutta Method")
plt.xlabel("Time in Seconds"); plt.ylabel("Theta Points in Radians")
plt.grid(True); plt.show()
