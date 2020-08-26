# DynamicNets

The code provides estimation and inference for *dynamic networks* measures introduced in the following papers

Baruník, J. and Ellington, M. (2020): *Dynamic Networks in Large Financial and Economic Systems*, manuscript [available here for download](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3651134) (July 2020) 

Baruník, J. and Ellington, M. (2020): *Dynamic Network Risk*, manuscript [available here for download](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3622200) (July 2020) 

## Time Horizon Dynamics

Note that current version of the codes works with 3 possible horizons of user's choice

Functions dynamic_networks.m and get_dynnet.m estimates time varying total network connectedness as well as directional connectedness, timing for a 4 variable system, 1832 time observations, and 100 simulations is around 230 seconds. This is for a Desktop PC with 64GB RAM  with 3.70 GHz 6-Core Intel Core i7 processor.

The toy data is daily and we provide an example of dynamic horizon specific network with horizons defined as
* short run: 1 - 5 days (up to one week)
* medium run: 5 - 20 days (week up to month)
* long run: 20 + days (more than month)
