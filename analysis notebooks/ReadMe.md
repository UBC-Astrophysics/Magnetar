Some analysis notebooks for IXPE data using the maximum likelihood techique outlined in https://arxiv.org/abs/2204.00140, and producing and using models with the Magnetar python package.

* [Reduce-4U.ipynb](Reduce-4U.ipynb)    : reduce the data for 4U 0142+61 including energy correction, period finding and folding
* [Reduce-CenX3.ipynb](Reduce-CenX3.ipynb) : reduce the data for Cen X-3 including source barycentring, period finding and folding
* [MCMC-Julia.jl](MCMC-Julia.jl) : Julia script to perform MCMC estimation of confidence regions on RVM parameters (uses photondata.txt produced by the reduce notebooks)
* [corner-plots-CenX3.ipynb](corner-plots-CenX3.ipynb) : Produce corner plots from the Julia MCMC code
