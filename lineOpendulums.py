import math
from math import pi

# Define some variables
lengths_meters = []  # Length in Meters, EAT THAT YOU IMPERIAL PEOPLE!
gamma = 60  # Our time in seconds
N = 51

# Where n is the index of the pendulum (starting at 0)
length = lambda n: 9.81 * (gamma/ (2*pi*(N + n)))**2

# 12 pendulums, indexed 0-11
for n in range(12):
    lengths_meters.append(length(n))

# Round our lengths to 3 decimals
lengths_meters_rounded = []
for length in lengths_meters:
    lengths_meters_rounded.append(round(length, 3))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# Set the backend of matplotlib to the 'inline' backend
#%matplotlib inline

# Set the pivot point
location = np.linspace(0, 0.056*11, 12)

# Create some functions
func_theta = lambda theta_i, n, t: theta_i * math.cos((2 * pi * t) / (60 / (51 + n)))
func_dispx = lambda l, theta: l * math.sin(theta)
func_dispy = lambda l, theta: -1 * l * math.cos(theta)

def init1():
    line.set_data([], [])
    return (line,)

def animate1(t):
    displacementArr = []
    for n in range(12):
        displacementArr.append(func_dispx(lengths_meters[n], func_theta(math.radians(10), n, t)))
    fig.suptitle("Top View at t = {:.3f} seconds".format(round(t, 3)), fontsize=12)
    line.set_data(location, np.array(displacementArr))
    return (line,)

# Set up the figure, the axis and the plot we want to animate
fig, ax = plt.subplots()
ax.set_xlim([0, location[11]])
ax.set_ylim([-0.06, 0.06])
line, = ax.plot([], [], linestyle='--', marker='o')
plt.xlabel("Position along the device, meters")
plt.ylabel("Horizonal displacement from the resting point (in meters)")
fig.suptitle("Top View at t = 0 seconds", fontsize=12)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.labelpad = 100

# Run the animation and save it as a gif with 60fps... 60 FPS MASTER RACE!
fps=60
anim = animation.FuncAnimation(fig, animate1, init_func=init1, frames=np.linspace(0, 60, 60*fps//2), blit=True)
anim.save('animation_top_bottom.gif', dpi=80, writer='imagemagick', fps=60)
plt.close()


fig, ax = plt.subplots()
ax.set_xlim([-0.15, 0.15])
ax.set_ylim([-0.4, 0])
line, = ax.plot([], [], linestyle='', marker='o')
plt.xlabel("Horizonal displacement from the resting point (in meters)")
plt.ylabel("Vertical displacement from the resting point (in meters)")
fig.suptitle("Top View at t = 0 seconds", fontsize=12)
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.xaxis.labelpad = 220
ax.yaxis.labelpad = 150

# Create some functions
def init2():
    line.set_data([], [])
    return (line,)

def animate2(t):
    x_pos = []
    y_pos = []
    for n in range(12):
        x_pos.append(func_dispx(lengths_meters[n], func_theta(math.radians(10), n, t)))
        y_pos.append(func_dispy(lengths_meters[n], func_theta(math.radians(10), n, t)))
    fig.suptitle("Head-On View at t = {:.3f} seconds".format(round(t, 3)), fontsize=12)
    line.set_data(np.array(x_pos), np.array(y_pos))
    return (line,)

# Run the animation and save it as a gif with 60fps... 60 FPS MASTER RACE!
fps = 60
anim = animation.FuncAnimation(fig, animate2, init_func=init2, frames=np.linspace(0, 60, 60*fps//2), blit=True)
anim.save('animation_front_rear.gif', dpi=80, writer='imagemagick', fps=60)
plt.close()
