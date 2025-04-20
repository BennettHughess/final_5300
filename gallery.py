# import necessary packages
import matplotlib.pyplot as plt
from matplotlib import animation as anim
from NBodyHandler import *
import numpy as np

##############################################################
# FIRST SIMULATION
# This first simulation will be a simple four-body system.

# initalize particles
mass_array = [1, 1, 1, 1]
position_array = [
    [-1, 0],
    [1, 0],
    [0, 1],
    [0, -1]
]
velocity_array = 0.5*np.array([
    [0, -1],
    [0, 1],
    [-1, 0],
    [1, 0]
])

# initialize the NBodyHandler
nbody = NBodyHandler(mass_array, position_array, velocity_array)

# create time arrays and solve using leapfrog4
tlist = np.arange(0, 17, 0.01)
t, vec = nbody.solve(tlist, method="Leapfrog4")

# plot with a mosaic : two plots on top, one on bottom
fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

# keep initialize the plot which gets updated
points_arr = []
pos_arr = []
for idx in range(0, len(mass_array)):

    if idx == len(mass_array)-1:
        label = r"$m_{test}$" + rf"$ = {mass_array[idx]}$"
    else:
        label = rf"$m_{idx} = {mass_array[idx]}$"

    point = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0],
                    '-', label=label)[0]
    points_arr.append(point)

    pos = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0], '.k', markersize=10)[0]
    pos_arr.append(pos)


# make the animation
def update(frame):
    for idx in range(0, len(mass_array)):
        points_arr[idx].set_xdata(vec[2*idx, :frame])
        points_arr[idx].set_ydata(vec[2*idx+1, :frame])
        pos_arr[idx].set_xdata([vec[2*idx, frame]])
        pos_arr[idx].set_ydata([vec[2*idx+1, frame]])

    return points_arr, pos_arr


ax.set_axis_off()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)

fig.tight_layout()

anim.FuncAnimation(fig=fig, func=update, frames=len(tlist),
                   repeat=True, blit=False
                   ).save(
    "images/four_body.mp4",
    writer="ffmpeg", fps=60
)

##############################################################
# SECOND SIMULATION
# Now let's do six bodies!

# initalize particles
mass_array = [1, 1, 1, 1, 0.1, 0.1]
position_array = [
    [-1, 0],
    [1, 0],
    [0, 1],
    [0, -1],
    [1, 1],
    [-1, -1]
]
velocity_array = 0.5*np.array([
    [0, -1],
    [0, 1],
    [-1, 0],
    [1, 0],
    [2, -2],
    [-2, 2]
])

# initialize the NBodyHandler
nbody = NBodyHandler(mass_array, position_array, velocity_array,
                     epsilon=0.025, G=1)  # add some buffer radius

# create time arrays and solve using leapfrog4
tlist = np.arange(0, 4, 0.01)
t, vec = nbody.solve(tlist, method="Leapfrog4")

fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

# keep initialize the plot which gets updated
points_arr = []
pos_arr = []
for idx in range(0, len(mass_array)):

    if idx == len(mass_array)-1:
        label = r"$m_{test}$" + rf"$ = {mass_array[idx]}$"
    else:
        label = rf"$m_{idx} = {mass_array[idx]}$"

    point = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0],
                    '-', label=label)[0]
    points_arr.append(point)

    pos = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0], '.k', markersize=10)[0]
    pos_arr.append(pos)


# make the animation
def update(frame):
    for idx in range(0, len(mass_array)):
        points_arr[idx].set_xdata(vec[2*idx, :frame])
        points_arr[idx].set_ydata(vec[2*idx+1, :frame])
        pos_arr[idx].set_xdata([vec[2*idx, frame]])
        pos_arr[idx].set_ydata([vec[2*idx+1, frame]])

    return points_arr, pos_arr


ax.set_axis_off()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)

fig.tight_layout()

anim.FuncAnimation(fig=fig, func=update, frames=len(tlist),
                   repeat=True, blit=False
                   ).save(
    "images/six_body.mp4",
    writer="ffmpeg", fps=60
)

##############################################################
# THIRD SIMULATION
# The remaining simulations are taken from the sick website
# https://observablehq.com/@rreusser/periodic-planar-three-body-orbits
# This one is the figure 8!

# initalize particles
mass_array = [1, 1, 1]
position_array = [
    [-1, 0],
    [1, 0],
    [0, 0]
]
velocity_array = np.array([
    [0.347113, 0.532727],
    [0.347113, 0.532727],
    [-0.694226, -1.065454]
])

# initialize the NBodyHandler
nbody = NBodyHandler(mass_array, position_array, velocity_array,
                     epsilon=0, G=1)

# create time arrays and solve
tlist = np.arange(0, 20, 0.05)
t, vec = nbody.solve(tlist, method="Leapfrog4")

fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

# keep initialize the plot which gets updated
points_arr = []
pos_arr = []
for idx in range(0, len(mass_array)):
    point = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0],
                    '-')[0]
    points_arr.append(point)

    pos = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0], '.k', markersize=10)[0]
    pos_arr.append(pos)


# make the animation
def update(frame):
    for idx in range(0, len(mass_array)):
        points_arr[idx].set_xdata(vec[2*idx, :frame])
        points_arr[idx].set_ydata(vec[2*idx+1, :frame])
        pos_arr[idx].set_xdata([vec[2*idx, frame]])
        pos_arr[idx].set_ydata([vec[2*idx+1, frame]])

    return points_arr, pos_arr


ax.set_axis_off()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)

fig.tight_layout()

anim.FuncAnimation(fig=fig, func=update, frames=len(tlist),
                   repeat=True, blit=False
                   ).save(
    "images/figure_8.mp4",
    writer="ffmpeg", fps=60
)

##############################################################
# FOURTH SIMULATION
# This simulation is "Brouke A1"

# initalize particles
mass_array = [1, 1, 1]
position_array = [
    [-0.9892620043, 0],
    [2.2096177241, 0],
    [-1.2203557197, 0]
]
velocity_array = np.array([
    [0, 1.9169244185],
    [0, 0.1910268738],
    [0, -2.1079512924]
])

# initialize the NBodyHandler
nbody = NBodyHandler(mass_array, position_array, velocity_array,
                     epsilon=0, G=1)

# create time arrays and solve
tlist = np.arange(0, 30, 0.01)
t, vec = nbody.solve(tlist, method="Leapfrog4")

fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

# keep initialize the plot which gets updated
points_arr = []
pos_arr = []
for idx in range(0, len(mass_array)):
    point = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0],
                    '-')[0]
    points_arr.append(point)

    pos = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0], '.k', markersize=10)[0]
    pos_arr.append(pos)


# make the animation
def update(frame):
    for idx in range(0, len(mass_array)):
        points_arr[idx].set_xdata(vec[2*idx, :frame])
        points_arr[idx].set_ydata(vec[2*idx+1, :frame])
        pos_arr[idx].set_xdata([vec[2*idx, frame]])
        pos_arr[idx].set_ydata([vec[2*idx+1, frame]])

    return points_arr, pos_arr


ax.set_axis_off()
ax.set_xlim(-3, 3)
ax.set_ylim(-2, 2)

fig.tight_layout()

anim.FuncAnimation(fig=fig, func=update, frames=len(tlist),
                   repeat=True, blit=False
                   ).save(
    "images/brouke_A1.mp4",
    writer="ffmpeg", fps=60*4
)

##############################################################
# FIFTH SIMULATION
# This simulation is "Brouke R9"

# initalize particles
mass_array = [1, 1, 1]
position_array = [
    [0.901558607, 0],
    [-0.6819108246, 0],
    [-0.2196477824, 0]
]
velocity_array = np.array([
    [0, 0.9840575737],
    [0, -1.6015183264],
    [0, 0.6174607527]
])

# initialize the NBodyHandler
nbody = NBodyHandler(mass_array, position_array, velocity_array,
                     epsilon=0, G=1)

# create time arrays and solve
tlist = np.arange(0, 30, 0.01)
t, vec = nbody.solve(tlist, method="Leapfrog4")

fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

# keep initialize the plot which gets updated
points_arr = []
pos_arr = []
for idx in range(0, len(mass_array)):
    point = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0],
                    '-')[0]
    points_arr.append(point)

    pos = ax.plot(vec[2*idx, 0], vec[2*idx+1, 0], '.k', markersize=10)[0]
    pos_arr.append(pos)


# make the animation
def update(frame):
    for idx in range(0, len(mass_array)):
        points_arr[idx].set_xdata(vec[2*idx, :frame])
        points_arr[idx].set_ydata(vec[2*idx+1, :frame])
        pos_arr[idx].set_xdata([vec[2*idx, frame]])
        pos_arr[idx].set_ydata([vec[2*idx+1, frame]])

    return points_arr, pos_arr


ax.set_axis_off()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)

fig.tight_layout()

anim.FuncAnimation(fig=fig, func=update, frames=len(tlist),
                   repeat=True, blit=False
                   ).save(
    "images/brouke_R9.mp4",
    writer="ffmpeg", fps=60*4
)
