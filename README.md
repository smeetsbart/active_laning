# active_laning
This repository contains the scripts for the simulations of self-propelled particle laning.
In this model, cells experience anisotropic friction and actively move in the direction of a polarization vector that undergoes
a specific type of diffusion.

Author: [Bart Smeets](https://github.com/smeetsbart)

Dependencies:
 - Python 3 and higher
 - Legacy [Mpacts](https://gitlab.kuleuven.be/mebios-particulate/mpacts) back-end, the particle-based simulation platform
 - Python numpy

The most important simulation script contained in this repository can be found in [python/simulation/anisotropic_field.py](python/simulation/anisotropic_field.py)

The script [python/simulation/anisotropic_field.py](python/simulation/anisotropic_field.py) is based on our paper on bidirectional laning from collective contact guidance,
which contains a combined experimental/theoretical/computational study on tissue laning as a result of anisotropic friction.
In this, we use the following equations of motion:

```math
\hat{\Lambda}_s\,\hat{\vec{v}}_i = \vec{p}_i - \sum_{j=0}^N \frac{\mathrm d \hat{U}_{ij}}{\mathrm d \delta_{ij}}
```
for the velocity $\vec{v}$, with $\hat{\Lambda_s}$ the adimensionalized substrate friction tensor, and $\hat{U}$ the adimensionalized potential strength, and
```math
\mathrm d \vec{p}_i = \left[ k_s \hat{\vec{v}}_i - (k_s+1) ||\hat{\vec{v}}_i||^2 \vec{p}_i  \right]\,D_p\,\mathrm d t + \sqrt{ 2 D_p }\,\mathrm d \vec{W}_i.
```
for the particle polarization $\vec{p}$. Here, $k_s$ is the self-reinforcement strength that parameterizes the velocity-polarity coupling and $D_p$ the diffusivity of the polarity.

```bash

$ cd scripts/
$ python3 anisotropic_field.py

```

You can use [Paraview](https://www.paraview.org/) for visualization of simulation output.





