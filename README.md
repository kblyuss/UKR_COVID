# Modelling the dynamics of COVID-19 in Ukraine

This folder contains **MATLAB** and data files needed to simulate the dynamics of COVID-19 in Ukraine.

Age-specific parameter values are contained in **AM3.dat** (1st column: *age*; 2nd column: % of *total cases*; 3rd column: % *mortality in each age group*; 4th column: % *cases requiring hospitalisation*; 5th colmumn: % *medics*; 6th column: % *cases returning from abroad*; 7th column: % *symptomatic cases*; 8th column: % *asymptomatic cases*).

Contact matrices for home, school, work, and other settings are contained in **ukr_home.txt**, **ukr_school.txt**, **ukr_work.txt**, and **ukr_other.txt**, respectively.

Simulation of an early stage of the epidemic, and a one-week forecast, compared to data: **short_term_forecast.m**

Simulation of a longer-term forecast, with possible increases of transmission by 10% or 20%: **long_term_forecast.m**

Simulation of interventions over a longer-term timescale: **quarantine_simulation.m**
