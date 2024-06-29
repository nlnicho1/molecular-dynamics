# Molecular Dynamics Simulation #


This simulation aims to demonstrate simplified interactions between noble gas atoms[^1].



## Description 

While interactions depend on molecular shapes and orientations, we are going to ignore this and assume spherical molecules. There is no exact formula to calculate the intermolecular forces, however, the *Lennard-Jones potential* gives us a close approximation: 

	V(r) = 4ε [(σ/r)^12 – (σ/r)^6]


![Graph_of_Lennard-Jones_potential](https://github.com/nlnicho1/molecular-dynamics/assets/80715072/0bd0ddd9-439d-4533-939b-32d88ec4bfc3)[^2]


The *Verlet Algorithm* was used to integrate Newton's Laws of motion:

  ![Verlet Algorithm](https://github.com/nlnicho1/molecular-dynamics/assets/80715072/0bd6e538-eced-4dcc-8d80-d7e5b2e2da6f)[^3]




[^1]: Simulation was inspired by Daniel V. Schroeder, PhD, [*Physics Simulations in Python*](https://physics.weber.edu/schroeder/scicomp/PythonManual.pdf)

[^2]: [Wikipedia - Lennard-Jones Potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential)

[^3]: [Wikipedia - Verlet Integration](https://en.wikipedia.org/wiki/Verlet_integration)
