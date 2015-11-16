README.TXT FOR COALESCENT SIMULATION

The file "AsexCoalescentSim.R" is the source code for the coalescent simulations used in the Hartfield, Wright and Agrawal 2015 paper "Coalescent times and patterns of genetic diversity in species with facultative sex: effects of gene conversion, population structure and heterogeneity". Simulations are written in R, which needs to be downloaded from  http://www.r-project.org/ before use.

EXECUTION

The script is meant to be run via a command-line. The command 'Rscript' therefore needs to be used. The general syntax is as follows: all commands are obligatory.

Rscript AsexCoalescentSim.R N GC Theta SwitchingType pLH pHL M Demes [Rates per deme] Reps

Each command is defined as follows:

- N is the population size of the ENTIRE population. So if there are d demes, there are N/d individuals per deme. Make sure that N/d is an integer!
- GC is the scaled rate of gene conversion, i.e. equal to 2Ng where g is the unscaled probability of gene conversion.
- Theta is the net neutral mutation rate, 4Nmu. Used to determine neutral diversity along the genealogy.
- 'SwitchingType' is a command to determine what kind of temporal heterogeneity in sex rates are present, if at all. The input has to be one of 0, 1, or 2.
	- 0 indicates transitions between two rates of sex over time ('Temporal change' in rates of sex). In this case, the initial rate of sex is determined based on a weighted mean of sex rates (so the more common transition is most likely to be the current state).
	- 1 indicates a stepwise change in sex (from a current state to ancestral state after a timepoint)
	- 2 indicates no change in sex rate over time ('constant-sex rate').
- pLH is one of the following, dependent on the value of 'SwitchingType':
	- If SwitchingType = 0, i.e. changing sex between two values, then pLH is the probability of switching from the low-sex state to high-sex state.
	- If SwitchingType = 1 above, i.e. stepwise change in sex, then this is the time in the past (in 2N generations) at which the stepwise change in sex arises.
	- If SwitchingType = 2 above, this parameter is not used.
- pHL is defined if SwitchingType = 0. It is the probability of switching from the high-sex state to low-sex state. Otherwise it is not used.
- M = 2Nm is the net migration rate between demes, for m between 0 and 1. If there only exists one deme, then M is internally set to zero as a precaution.
- Demes are the number of subpopulations present. Set to 1 to define a single panmictic population.

[Rates per deme] are a set of 4*Demes parameters, which together determine how sex is defined within each deme. Together they allow the user to define both temporal and spatial changes in sex rates over time, if desired. For each deme four parameters are defined as follows.
- Number of paired, within-individual samples present in that deme at the present time (note that each sample equates to two separate samples, one from each gene copy within an diploid individual)
- Number of unique, single samples from individuals in that deme at the present time
- The 'Low' rate of sex in that deme
- The 'High' rate of sex in that deme

If SwitchingType = 1 (stepwise change) then the 'Low Sex' parameter is the rate of sex in the current time, while "High Rate" is the ancestral rate of sex following a stepwise change in rates.
If SwitchingType = 2 (no temporal heterogeneity) then only the 'Low' rate is used.

Finally, 'Reps' is the number of time to repeat the simulation, to produce 'Reps' genealogies.

OUTPUT

Output is very similar to that produced by the simulation program MS, to allow compatibility between the two. The following files are produced:
- "Seed.out" containing the random seed used at the start of the simulation
- "Trees.out" with NEWICK printout of the trees produced for each repetition of the simulation
- â"CoalTimes.dat", a table listing the occurrence of the n-1 coalescent times for each repetition of the simulation (for n the number of samples present)
- A folder "Mutations", where each file lists the mutation spectrum of the population. Each row represents an individual mutation, with position between 0 and 1 (column 1), followed by whether it is absent (0) or present (1) in each of the n samples.

Information is also piped to STDOUT at the end of the simulation, so it is wise to redirect this to another file during execution.

EXAMPLES

The simplest case is a single population with a single rate of sex. The command for such a case would be similar to this:

Rscript AsexCoalescentSim.R 10000 0.001 5 2 0 0 0 1 20 10 0.01 0 1000

For a diploid population of size 10000; a gene conversion rate of 0.001; mutation rate of 5; "SwitchingType" = 2 indicating no temporal heterogeneity in sex rates, so the next two parameters are 0; 0 migration since there is only 1 deme; 20 paired samples and 10 single samples within the population; rate of sex = 0.01; simulation repeated 1000 times.

If we have two demes with different rates of sex within each, say 0.001 in one and 1 in the other and 10 paired samples in each, then after setting d = 2 there are now eight parameters to define. The following setup also assumes a migration rate M = 2Nm = 0.1.

Rscript AsexCoalescentSim.R 10000 0.001 5 2 0 0 0.1 2 10 0 0.001 0 10 0 1 0 1000

Now say we want to investigate temporal heterogeneity in sex within one population. The fourth parameter ("SwitchingType") is set to zero to denote this. The fifth and sixth parameters are set to pLH and pHL, which are 0.001 and 0.01 in this case (which are rather high values to use!). The population contains 25 paired samples, switching from low-sex rate of 0.001 to high-sex rate of 1.

Rscript AsexCoalescentSim.R 10000 0.001 5 0 0.001 0.01 0 1 25 0 0.001 1 1000

Below the commands are altered to consider a stepwise change in sex rates: SwitchingType is set to 1, and the fifth parameter equals 2, meaning that the transition will take place 4N generations in the past. The current rate is 0.001, switching to 1.

Rscript AsexCoalescentSim.R 10000 0.001 5 1 2 0 0 1 25 0 0.001 1 1000

Finally, a complicated example. Say we wish to investigate a three-deme population, where the rates of sex are different in each deme AND change over time. Let the transition probability over time be equal to 0.0001 both ways (i.e. every 10000 = N generations), M = 1. In deme 1 the low rate of sex equals 0.001 and the high rate is 1. In deme 2, the low rate is 0.5 and the high rate is 1. Finally in deme three, the low rate is 0.00001 and the high rate is 0.0001. There are 10 paired samples in each deme. The command line to set up such a simulation is as follows:

Rscript AsexCoalescentSim.R 10000 0.001 5 0 0.0001 0.0001 1 3 10 0 0.001 1 10 0 0.5 1 10 0 0.00001 0.0001 1000

Note that when a rate of sex changes over time, it does so in all demes simultaneously.

LOW SEX SIMULATION

Also provided is a variant of the simulation, 'AsexCoalescentSim_LowSex.R', that assumes all sex is low (that is, N*sigma is of order 1). For this simulation the population size does not have to be defined, as in traditional coalescent simulations. This program is also run via the command line using the 'RScript' command:

Rscript AsexCoalescentSim_LowSex.R GC Theta SwitchingType pLH pHL M Demes [Rates per deme] Reps

Each command is defined as follows: highlighted are changes in parameters compared to the full simulation. The main change to bear in mind is that parameters (rather than switches) are all scaled to 2N or 4N, which is not the case with the complete simulation.

- GC is the scaled rate of gene conversion, i.e. equal to 2Ng where g is the unscaled probability of gene conversion.
- Theta is the net neutral mutation rate, 4 N mu. Used to determine neutral diversity along the genealogy.
- 'SwitchingType' is a command to determine what kind of temporal heterogeneity in sex rates are present, if at all. The input has to be one of 0, 1, or 2. The definitions are the same as in the full simulation.
- If SwitchingType = 0, i.e. changing sex between two values, then pLH is the *scaled rate (in time units of 2N generations)* of switching from the low-sex state to high-sex state. *For example, a value of pLH = 2 will mean that the rate of sex will change every N generations, on average. This is different from the previous simulation, where this value was defined as a probability instead*. The other cases of pLH (for SwitchingType = 1, 2) are defined as in the full simulation.
- pHL is defined if SwitchingType = 0. It is the *rate (in 2N generations)* of switching from the high-sex state to low-sex state. Otherwise it is not used.
- M = 2Nm is the net migration rate between demes, for m between 0 and 1. Defined as in the full simulation.
- Demes are the number of subpopulations present. Set to 1 to define a single panmictic population.

[Rates per deme] are a set of 4*Demes parameters, which together determine how sex is defined within each deme. *Here, the rates of sex are defined in units of 2N generations. Hence a value of 2 equates to an unscaled rate of sex of 1/N in the previous simulations.*

A simple example:

Rscript AsexCoalescentSim_LowSex.R 0 5 2 0 0 0 1 25 0 2 0 1000

This command runs 1000 coalescent simulations with 25 within-indivdual initial samples, with 2Nsigma = 2 (i.e. unscaled sigma = 1/N), no gene conversion and theta = 5.

CONTACT

Comments should be sent to Matthew Hartfield (matthew.hartfield@utoronto.ca).
