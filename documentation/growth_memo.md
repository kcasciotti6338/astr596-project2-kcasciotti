# Growth Memo
## What you built and why

I built an N-body simulator. Within this, I have a class called N_body, which has mass, positions, velocities, N, dimensions, radial densities, softening factor, radius, kinetic energy, potential energy, total energy, virial ratio, and acceleration as attributes. I started out with just mass, position, velocity, and acceleration, and added attributes as I needed them for implementation. Within here, is the acceleration and energy functions, since those are specific to a certain instance for a system, and can't be calculated as a whole after integration. 

I then have my integration methods module, which contains a large function called integrate, which contains and calls the other integration methods. This can also save only certain snapshots of the integration to save on memory/time. The module also contains a dataclass of results, which holds the results of the integration, like times, positions, velocities, energies, and other plotting parameters. I used a dataclass for this because I wanted easy access to my results for plotting, without adding any heavy get functions to my N_body class, which would slow things down. I did some googling and found that since I only wanted to use this after integration, there was no need to make increasingly complex numpy arrays, or dictionaries or another class, which would lead to more complexity and repetitiveness and open to bugs. This results dataclass has not caused me a single problem, it did exactly what it needed to without taking up much space, both in the code and in my mind. 

I also have a sampling methods module, which contains everything needed to initialize a system. There are sun-earth and sun-earth-jupiter functions that just hold those initial conditions that aren't changing between tests. The main thing in the module is the NBody sampling function. This takes N, a and spits out a fully initialized N-body object. Within here, there are the Kroupa IMF and Plummer sphere samplers. Originally, these were separate functions, but once I had them working, I had no intention of ever using them outside of the N-body sampler, so I put them inside. 

I have my plotting module, which is pretty straight forward. I also made sure to save all of my old testing scripts, so there's a separate folder for old tests and bad plots. 

At the end, I made a final analysis module with all of the required plots and tests detailed in functions, which I can call individually later as I compile my best final plots. 

## Challenges you faced (everyone has them!)

I struggled a lot with the vectorization. I spent probably 6-8 hrs just writing out pseudo code, reviewing matrix operations, and actually attempting to implement. I had no trouble conceptualizing the loops, it was just figuring out what could broadcast to where and which axes I need to sum or multiply over. I'm glad I did this first, though. I often get stuck trying to make the first step perfect and then get frustrated when I have to delete half of it because it has no place. This way, I forced myself to start with the big picture and work my way in. 

I also spent way too long trying to get my plots to do what I wanted them to do. My plotting module ended up a mess. I kept plotting and wanting better, so I'd look up the matplotlib examples, and try to implement something from like 10 different examples without changing my base structure first, so then data types and variable names kept clashing and things kept erroring out and not being any better than how I started. Eventually, I went back through and tried to cut out everything that was unnecessary, but it's still a mess. 

## How you used and verified AI assistance

There were a few times where things were going wrong, and I could not find the bug after like 30-45 mins, so I sent my code over to ChatGPT, and it was able to find my mistakes. The most prominent ones were when I forgot to not use G in cgs units and when I accidentally added instead of subtracted in one of the plummer sphere equations. 

I also tried to use it to clean up my plotting module, but I didn't think I was allowed to copy and paste from AI, so I just asked it where my variables were getting tangled up and how to merge ideas from various examples, and it helped a bit. There were a couple of things it suggested that seemed like weird bandaid fixes, like adding a random empty class into the middle of my plotting function, but I couldn't figure out another solution, and that worked, so it's there. I'll definitely start my plotting module from scratch next time and go in slowly with a plan, instead of chaotically adding things left and right. 

## What excited or surprised you

I was surprised at how finicky every parameter was to the results. I messed around with a bunch of different radii, total time, time steps, and N values, and my outputs looked completely wrong a lot of the time, simply because I had the radius too big for that sample, or my step size was too big, or too small, so it took way too long to compute. 

## What you’d explore next

I would like to dive deeper into optimization. I took an optimization class in C++ last year, so I know the basics, but knowing theory and actually putting it into practice are different things. I always forget to think about it until the end, when my code takes forever to run, so I think I'll force myself to include it in my pseudo code next time, so that I don't get stuck with computationally heavy structures.
