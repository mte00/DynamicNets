# DynamicNets

The code has been developed in Julia as a code accompanying the Barunik and Ellington (2023) and provides an estimation and inference for *dynamic networks* measures introduced in

Baruník, J. and Ellington, M. (2023): *Persistence in Financial Connectedness and Systemic Risk*, European Journal of Operational Research, forthcoming, manuscript [available here for download](https://ideas.repec.org/p/arx/papers/2007.07842.html) (revised Sep 2023)

## Time Horizon Dynamics

Note that current version of the codes works with 3 possible horizons of the user's choice

Functions dynamic_networks.m and get_dynnet.m estimates time varying total network connectedness as well as directional connectedness, timing for a 4 variable system, 1832 time observations, and 100 simulations is around 230 seconds. This is for a Desktop PC with 64GB RAM  with 3.70 GHz 6-Core Intel Core i7 processor.

The toy data is daily and we provide an example of dynamic horizon specific network with horizons defined as
* short run: 1 - 5 days (up to one week)
* medium run: 5 - 20 days (week up to month)
* long run: 20 + days (more than month)

## Main Files

Dynamic_Nets_Master.m is the master file
dynamic_networks.m is the function with the TVP VAR estimation. See inputs and outputs in the file.  
Within this function, you may want to change the variables: 
* w which denotes the kernel width (default is set to 8).
* HO which denotes the horizon you compute the wold decomposition for (default set to 10+1 for speed, for applications you should set to large value such as 100. This will cause computation time to increase).

get_dynnet.m is the function that computes the time-frequency network measures. See inputs and outputs in the file.
Within this function, you may want to change the variables:
* d1 which determines long-term definition (default is set to >20-days)
* d2 which determines medium-term definition (default is set to 5-days to 20-days)
* d3 which determines short-term definition (default is set to 1-day to 5-days)
### Note that d1, d2, d3 are dependent on the frequency you observe your data and will depend on the application at hand. You will need to look at omeg2 to determine your horizons.

### Update 13/05/2021 Code now provides ADJACENCY MATRICES across horizons and aggregated over horizons. From this you can readily compute pairwise connections. 
function get_dynnet.m provides the adjacency matrices at each horizon and also aggregated across all horizons (Diebold Yilmaz, 2014).
Dynamic_Nets_Master.m is the master file and provides the user with the NxNxT adjacency matrices which is the posterior median over draws from the QBLL methodology. 
dynamic_networks.m updated to save adjacency matrices at each draw from posterior distribution.
You will need to download these scripts again to estimate.
