# swirlab

swirlab is a Matlab code for modeling shortwave infrared (SWIR) 
radiative transfer and for retrieving atmospheric profiles 
from FTIR data.

swirlab needs pre-computed absorption cross sections to work.
Some example cross sections for CH4 can be downloaded from: https://goo.gl/442NoI

Another option is to compute them from HITRAN parameters and Q-values: https://goo.gl/5TbkzX

swirlab uses MCMC Matlab toolbox by Marko Laine. 
It is available at: https://mjlaine.github.io/mcmcstat/

test_ch4_retrieval.m is a simple example how to run swirlab. 

Use matlab version 2016b or newer!
