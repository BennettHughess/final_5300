# Physics 5300 final: N-body simulation
**Please read all of this before diving into the project!**

This is my final project for Physics 5300, where I've implemented and discussed an n-body simulation. This readme will explain the contents of the project, an overview of how I've implemented the rubric, and installation instructions.

## Getting started

1. To begin, either clone this repository or download it as a zip file.
```terminal
git clone https://github.com/BennettHughess/final_5300
```

2. Install the required Python packages, included in `requirements.txt`:
```terminal
pip install -r requirements.txt
```

3. Install [FFmpeg](https://ffmpeg.org/). If you are on MacOS and have Homebrew installed, `brew install ffmpeg` works. **If you don't install ffmpeg, matplotlib will not be able to save any animations.** This means that your computer won't be able to run the plotting code. 
    - However, all of the animation videos have been pregenerated, so you can skip this step if you just want to view the `.ipynb` file -- it'll automatically load the videos and images from the `images` directory.

4. Open `nbody.ipynb` using whatever method you want (VSCode, `jupyter notebook`) and execute the code! **Note: it may take several minutes for the code to execute, because FFmpeg has to do a good deal of work to save all of the animations.**

## Contents

This repository contains the following files:

- `nbody.ipynb` is the main file, and is an explanatory notebook describing the problem and implementing all of the plotting code.
- `NBodyHandler.py` contains the class `NBodyHandler`, and is used to perform the physical calculations.
- `/images` is the directory that contains all of the media files in the project, including videos. `nbody.ipynb` generates almost all of the files in this directory, but the directory is prepopulated from when I ran the code.
- `nbodytest.ipynb` was used for testing `NBodyHandler.py`, and can be safely ignored.
- `gallery.py` is used for creating the media in the gallery in this readme.

## Rubric

The final is out of 300 points -- I'll explain where I believe I've met the requirements on the rubric.


1. (Required, 50pts) Physical Problem
    - `nbody.ipynb` contains a diagram of the setup.
    - It walks through the derivation of the physical problem and calculates the equations of motion.
    - There is also a section dedicated to the symmetries and assumptions of the problem.
2. (Required, 50pts) Visualizations
    - `nbody.ipynb` contains many different plots of the motion of the bodies for different kinds of motion. It also contains helpful plots describing other features I wanted to emphasize.
    - I added a section containing vector field plots for this system as well.
    - Finally, I added a "gallery" section to this readme featuring other plots I made.
3. (Required, 25 pts) Submission of project code on github
    - That's this! :)
4. Exploration of Symmetries and Conserved quantities (50pts)
    - I have a section dedicated to the conservation of energy and momenta.
    - I also have a separate section dedicated to the chaos exhibited in the system.
5. Animation (75pts)
    - I used `matplotlib`'s `FuncAnimation` to make many different animated plots in this project, which are additionally saved in the `/images` directory.
6. Utilization of advanced numerical approaches (50 - 75pts)
    - I implemented the Euler and leapfrog integrators in the code and explored their differences. I also implemented a custom fourth-order leapfrog integrator using Yoshida coefficients, and I compared it to an adaptive RK45 algorithm.
    - I explored the advanced and disadvantages of adaptive vs symplectic integrators for this kind of problem.
7. Markdown Documentation of the code: (25pts)
    - `nbody.ipynb` has a good amount of markdown which explains the scenarios I'm simulating and adds commentary to the different concepts I'm exploring.

## Possible issues

- If `nbody.ipynb` is not running, make sure that all of the required dependencies are satisfied. 
    - Also make sure that the Jupyter kernel is restarted -- it shouldn't matter if the notebook is run sequentially, but it is good practice.
- If `matplotlib` is failing to make animations, make sure that FFmpeg is installed (optionally, you can rewrite the code to use `Pillow` to write the animations into `.gif`s with a little work).
    - If the animations still aren't working, you can always just look at the images and videos already contained in this repository which have been pregenerated on my machine.

## Gallery
These are just some pretty animations I generated.

[video](https://github.com/user-attachments/assets/a6a344ab-ef6a-4c16-8c5a-72c00aac862f)

[video](https://github.com/user-attachments/assets/8517aa46-9cdf-4466-afc2-ae100b4bbd83)

These next animations are initialized with conditions from [this helpful website](https://observablehq.com/@rreusser/periodic-planar-three-body-orbits).

[video](https://github.com/user-attachments/assets/21a2ac3d-2a84-45d1-9d00-9052a406945d)

[video](https://github.com/user-attachments/assets/edb82690-16ec-41b4-aa16-e47458d9f571)

[video](https://github.com/user-attachments/assets/858e9cbc-93f4-442c-b0d7-fb2862b3442f)





