# MATLAB code for modelling root water uptake in heterogeneous soil

This script includes model code of the following article:

*"Spatial heterogeneity enables higher root water uptake in dry soil but protracts water stress after transpiration decline: a numerical study"*

by Patrick v. Jeetze, Mohsen Zarebanadkouki, Andrea Carminati

## Descriptions

This script loads a predefined geometry given here as "SmallHexW51.png". This geometry contains polygons with different
gray values. These polygons represent regions of soil with varying hydraulic properties. In the first step the gray
value of these pre-drawn polygons is randomly distributed across the flow domain (using a script called 'CreateRdmDomain').
To do so we consider only ten randomly distributed regions. Then in the next step a randomly distributed hydraulic properties
is assigned to each region. The randomly distributed hydraulic properties(for only ten regions) are generated using an R script
(place in folder called 'KhEnsemble'). These steps reconstruct a flow domain will ten regions randomly distributed with randomly
redistributed hydraulic properties. Then the flow equation (Richards equation without gravity term) is solved by imposing a known
matric potential as initial condition, no flux on top and bottom, fixed matric potential on the right side (representing the bulk soil)
and a variable flux on the left side (representing the root surface). Flux at the root surface increases stepwise till soil water
potential at the root surface reach a critical value (set here ca. -10000 cm) and then it decreases to zero.


## To Run the Script:

The provided scripts can be used as follows:
open the script called MainScript.m and run it.
