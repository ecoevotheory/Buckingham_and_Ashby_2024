This ReadMe file contains information on files associated with the paper: "Coevolution of age-structured tolerance and virulence" by Buckingham & Ashby.

All files are provided "as is", without any express or implied warranty. The authors accept no liability for damages of any kind. 

Author: Lydia Buckingham, University of Bath
Contact: ljb74@bath.ac.uk
Date: 07/11/23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMMENTS

As a code naming convention, "LL" refers to the model scenario where lifelong tolerance evolves and "JL" refers to the scenario where only juvenile tolerance evolves (with adult tolerance held fixed).

MEX files written in C# must be compiled before use, using " mex codename.c ".

Be aware that some of the code relies on numerical approximations and so results may not be exact. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOURCE CODE FOR FIGURES

Figure 1:   "tradeoff_plots_fig1.m"		- Source code for plotting example trade-off functions												- written in matlab (R2019b). 

Figure 2:   "effect_of_lifespan_fig2.m"		- Source code for plotting the effect of lifespan on tolerance and virulence coevolution							- written in matlab (R2019b). 

Figure 3:   "cycling_simulation_fig3.m"		- Source code for plotting simulations of cycling in tolerance and virulence traits								- written in matlab (R2019b). 

Figure 4:   "cycling_curve_fig4.m"		- Source code for plotting the effect of lifespan on cycling											- written in matlab (R2019b). 

Figure S1:  "lifespan_mutationrate_figS1.m"	- Source code for determining the effect of lifespan and mutation rates on the incidence of bistability and cycling				- written in matlab (R2019b). 

Figure S2:  "lifespan_costs_figS2.m"		- Source code for determining the effect of lifespan and trade-off costs on the incidence of bistability and the evolution of full tolerance	- written in matlab (R2019b). 

Figure S3A: "effect_of_lifespan_figS3A.m"	- Source code for plotting the effect of lifespan with baseline transmissibility held fixed (and disease prevalence varying)			- written in matlab (R2019b). 

Figure S3B: "effect_of_lifespan_figS3B.m"	- Source code for plotting the effect of lifespan with disease prevalence held fixed (and baseline transmissibiltiy varying)			- written in matlab (R2019b). 

Figure S4:  "phaseplanes_figS4.m"		- Source code for drawing phase planes for different host and parasite relative mutation rates							- written in matlab (R2019b). 

Figure S5:  "cycling_curve_phaseplanes_figS5.m"	- Source code for plotting the effect of lifespan on cycling, illustrated by phase planes							- written in matlab (R2019b). 

Figure S6:  "hopf_bifurcation_figS6.m"		- Source code for drawing an Argand diagram of eigenvaluesof the full coevolutionary system							- written in matlab (R2019b). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTIONS USED IN SOURCE CODE

"LL_singstrat_function.m"			- Function which finds co-singular strategies and classifies their stability (in the lifelong tolerance scenario)				- written in matlab (R2019b).
"JL_singstrat_function.m"			- Function which finds co-singular strategies and classifies theit stability (in the juvenile tolerance scenario)				- written in matlab (R2019b).
"LL_simulation_function.m"			- Function which runs an evolutionary simulation for a fixed number of timesteps (in the lifelong tolerance scenario)				- written in matlab (R2019b).
"JL_simulation_function.m"			- Function which runs an evolutionary simulation for a fixed number of timesteps (in the juvenile tolerance scenario)				- written in matlab (R2019b).
"LL_fitness_gradients.m"			- Function for calculating the fitness gradients at different values of the coevolving traits (in the lifelong tolerance scenario)		- written in matlab (R2019b).
"JL_fitness_gradients.m"			- Function for calculating the fitness gradients at different values of the coevolving traits (in the juvenile tolerance scenario)		- written in matlab (R2019b).
"LL_classification_function.m"			- Function which determines the evolutionary and convergence stability of a co-singular strategy (in the lifelong tolerance scenario)		- written in matlab (R2019b).
"JL_classification_function.m"			- Function which determines the evolutionary and convergence stability of a co-singular strategy (in the juvenile tolerance scenario)		- written in matlab (R2019b).
"LL_find_parasite_singstrat_derivative.m"	- Function which finds the derivative of the parasite singular strategy with respect to the host trait (in the lifelong tolerance scenario)	- written in matlab (R2019b).
"JL_find_parasite_singstrat_derivative.m"	- Function which finds the derivative of the parasite singular strategy with respect to the host trait (in the juvenile tolerance scenario)	- written in matlab (R2019b).
"JL_sim_traj_function.m"			- Function for simulating an evolutionary trajectory given initial conditions (in the juvenile tolerance scenario)				- written in matlab (R2019b).
"fitgrad_signchange_function.m"			- Function which determines when the fitness gradient changes sign										- written in matlab (R2019b).
"singstrats_at_0or1.m"				- Function which determines cosingular strategies where one of the coevolving traits takes its maximum or minimum value				- written in matlab (R2019b).
"endemic_equilibrium_function.m"		- Function for determining the endemic equilibrium of the ecological system									- written in matlab (R2019b).
"LL_eco_dynamics_function.c"			- Function which runs ecological dynamics for a single evolutionary timestep (in the lifelong tolerance scenario)				- written as a MEX file in C#. 	
"JL_eco_dynamics_function.c"			- Function which runs ecological dynamics for a single evolutionary timestep (in the juvenile tolerance scenario)				- written as a MEX file in C#. 	
"ode_function.c"				- Function for running a system of ODEs to find its equilibrium											- written as a MEX file in C#. 					 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OTHER CODE FILES

"LL_fitgrad_expressions.m"			- Code for determining analytical expressions for the derivatives of the invasion fitness (in the lifelong tolerance scenario)			- written in matlab (R2019b). 
"JL_fitgrad_expressions.m"			- Code for determining analytical expressions for the derivatives of the invasion fitness (in the juvenile tolerance scenario)			- written in matlab (R2019b). 

See code for full description and instructions for use. 