clear all; load Model_setup3;

option = {'1718','uncon'};

Get_initial_guesses;
Do_MCMC2;
Get_MCMC_fitresults;