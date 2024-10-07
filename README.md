# QuantumSim

QuantumSim is a personal project targetting the need to simulate non-relativistic quantum systems.in 
The philosophy behind the code is that it must work in real-time, with an almost instantaneous setup of the initial wave-function and an optimized evolution algorithm.
Instead of targetting the most accurate description of reality, it targets the possibility to use those simulations in a fiction environment like in a video game or graphics animations.
Thus, the default numerical time evolvers are lightweight and stable in all configurations, but one can chose to use its own evolver if it requires more precision.

The purpose is to numerically solve the Schrödinger's equation:

$$\hat{\mathcal{H}} \psi = i \hbar \frac{\partial \psi}{\partial t}$$

The user sets the geometry of the system, included parameters and boundaries, and the potential at each timestep.
Than the simulation starts or can be sinchronized to an external clock, like an animation fps.

The project is divided into sections, each referring to the physical dimension they target: one-dimensional, bi-dimensional or three-dimensional.
Each dimension requires different solutions in order to obtain a stable simulation.

## Installation

Coming soon ...

## Usage

### One-dimensional
Coming soon ...

### Two-dimensional

## Applications
An example of application can be found here: [https://clone4004.itch.io/out-of-the-box](https://clone4004.itch.io/out-of-the-box)
The video game is incomplete if you consider its playability, but the wave the cat can spread is entirely simulated with the one-dimensional *quantumsim* framework.
