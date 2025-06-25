# Progression-of-retinal-detachment

A MATLAB code for simulating a fluid-structure interaction (FSI) model to study the progression of retinal detachment.
A detailed description of the model formulation and discretization can be found in the supplementary paper.

Main Files
MainFlow.m – This is the main script that runs the simulation and performs the time iterations.

OutPut.m – Defines all model parameters and calls the MainFlow function.


Description of Other Function Files
LagStructure.m – Computes the initial geometry of the retina.

AdhesionForce.m – Computes the adhesion force (F_h) binding the neural layer (NL) to the retinal pigmented epithelium (RPE).

BendForce.m – Computes the bending force (F_b) in the model.

BondDispW.m – Computes the displacement of the retina from its resting position [i.e., w = (X-X^e)].

ClampForce.m – Computes the total force at the clamped edge and the CFL condition to monitor code stability.

Delta.m – Computes the Dirac delta function linking the fluid domain to the retinal domain.

DetachRate.m – Computes the detached length to monitor the progression of detachment.

Falpha0.m – Computes the first term of the fluid force (F_l) on the retina using Goldstein’s feedback law.

FluidForce.m – Computes the force the retina exerts on the fluid.

InitRhib.m – Computes the initial distribution of the adhesion protein density (ρ_b).

InterpPoly.m – Performs interpolation for points outside the retinal domain needed for bending force computation.

LagForce.m – Computes the Lagrangian force (F_l) exerted by the fluid on the retina using Goldstein’s law.

P1FirstIter.m – Solves the Navier–Stokes equations explicitly for the second initial condition in the P₂ projection method.

P1RHS.m – Computes the right-hand side vector used in P1FirstIter.m

P2RHS.m – Computes the right-hand side vector needed for the P₂ projection method.

PLapMat.m – Computes the Laplacian matrix used to solve the Poisson equation for pressure correction.

RhoB.m – Computes the adhesion protein density (ρ_b) at each time step.

Tension.m – Computes the tension force in the retina.

UV_LapMat.m – Computes the Laplacian matrices for the u and v components of the fluid velocity.

Uib.m – Computes the fluid velocity at the immersed boundary.

UpdateLag.m – Updates the position of the retina at each time step.
