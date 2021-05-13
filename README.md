# DynamicNets

The code provides estimation and inference for *dynamic networks* measures introduced in the following papers

Baruník, J. and Ellington, M. (2021): *Dynamic Networks in Large Financial and Economic Systems*, manuscript [available here for download](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3651134) (February 2021) 

Baruník, J. and Ellington, M. (2020): *Dynamic Network Risk*, manuscript [available here for download](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3622200) (July 2020) 

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
