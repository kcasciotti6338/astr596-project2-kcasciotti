# Project Description: 

This project was to simulate a gravitationally bound N-body system. I also implemented four different iterative integrators (Euler's, Runge-Kutta 2 and 4, and Leapfrog) to demonstrate how a seemingly fine numerical method can become chaotic when scaled larger. 

To test this system, I started small with 2 bodies in 2 dimensions, the Earth orbiting the Sun. Then, I added Jupiter, and finally added 3 dimensions and 10+ bodies. 

# Installation instructions:

ipython

packages: numpy, matplotlib, imageio, Pillow

## Plots and Graphs:

I made a standard diagnostic that works for all systems. On a 2x2 figure, it shows the orbit or initial snapshot, the energy components over time, the relative energy error over time, and the virial ratio over time. 

I also included a histogram of the Kroupa IMF, snapshots of large systems, a log-log time vs N analysis plot, and plots showcasing very large N values. 

# Usage examples:

body_2(Tt=10, dt=0.01, snaps=200, gif=True)

This runs the Earth-Sun system for 10 years, with a stepsize of 0.01 yrs, and saves 200 evenly spaced snapshots throughout the total time. It goes through each of the four integration methods, plotting euler's on its own standard diagnostic figure and the other three on the same diagnostic figure for meaningful comparisons. Then, a gif of the orbit over one year will be saved as well, for each integration method. 

Each N-body demonstration has a similar style of analysis. 

# Key results summary:

For the integrator behavior, each method acted generally as expected. Euler's was never accurate enough to make even half and orbit before drifting into an uncontrolled spiral. RK2 was better but showed drift even at small scales. RK4 was much better, it showed accurate orbits with a very slow drift that was unnoticeable until the number of bodies and number of steps got larger. Leapfrog, as a symplectic model, did the best. It was very stable in terms of energy, with some drift when I started to push into much higher body and step counts. However, this was mostly due to other factors, such as the ratios between step size, softening, the number of bodies, and radius. 

# Acknowledgment of collaborators:

Throughout the project, I collaborated with my classmates to help conceptualize things and bounce ideas around. We talked the most about what to vectorize and what to leave alone, how the Plummer sphere works with re-centering, and ways to improve plots. We all share offices and classes, so there was a constant stream of communication and collaboration, so I couldn't name anyone or anything too specific. We never shared physical code, but did share sections of pseudo code and wrote out equations together. 

I also went to office hours, where a parenthesis bug was found in my force calculations, and got advice on improving my plots. 

I used NumPy, Matplotlib, Imageio, and Pillow.

In terms of internet and AI usage, I used ChatGPT to help debug, I used a lot of matplotlib examples to build my plotting functions, and I used Grammarly to spell check my written reports.
