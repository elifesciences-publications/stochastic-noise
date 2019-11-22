File descriptions:

COMMON PARAMETER AND FUNCTION DEFINITION FILES

ParameterFile.m - Set of default parameters used to run gene expression simulations.

FanoFactor.m - Function to calculate the Fano factor given single cell protein counts for output from 2 alleles

---------------------------------------------------------------------------
FILES REQUIRED TO SIMULATE A PAIR OF INDEPENDENT ALLELES

Simulation_Driver.m - Driver file to run simulations for independent alleles. Contains initial conditions for simulations such as number of cells and total simulation time. Also sets the range for rate K-on which determines the total range of protein output. Calls on TwoStatePromoter.m and ParameterFile.m to execute.

TwoStatePromoter.m - Function that implements Gillespie's Stochastic Simulation Algorithm for a single allele, given a set of rate parameters. Outputs stochastic mRNA and protein counts over time.

Analysis_code.m - Extracts allele-wise protein outputs for simulated cells at desired time-points. Sorts cells into bins according to total protein count and calls on FanoFactor.m to calculate bin-wise Fano factors. Plots single cell protein count correlation between alleles as well as Fano Factor vs protein for the simulation.

---------------------------------------------------------------------------
FILES REQUIRED TO SIMULATE A PAIR OF INTERACTING ALLELES

Simulation_Driver_interacting.m - Driver file to run simulations for interacting alleles. Calls on TwoStatePromoter_interacting.m and ParameterFile.m to execute.

TwoStatePromoter_interacting.m - Function that that implements Gillespie's Stochastic Simulation Algorithm for a single cells's stochastic mRNA and protein output for two interacting alleles.

Analysis_interacting.m - Data analysis code for interacting alleles. Plots single cell protein count correlation between alleles as well as Fano Factor vs protein for the simulation.


