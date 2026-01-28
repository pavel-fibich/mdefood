# mdefood
Simulation code and data, observed data for MDE foodwebs manuscript : "Does the mid-domain effect shape interaction networks along
environmental gradients".

Simulations:
* To run simulations, run simfood2.R,
* Outputs of the simulations (random and decreasing specialization) are in simout/,
* To reproduce Fig. 3,4,S3,S4 run figSims.R,
* To reproduce Fig. S1 run figs1sche.R.

Real food webs:
* Real food webs data (Myrmecophytic ant-plant dataset, Temperate plant-pollinator dataset, Tropical wet-season plant-pollinator dataset, and Tropical dry-season plant-pollinator dataset) are in realfood/,
* To get simulated (shuffled) data from real datasets run realfood.R,
* Outputs of shuffled data are in realfood/netsgrad1000noi_* files,
* To reproduce Fig. 5, S6 of real food webs run figReals.R.