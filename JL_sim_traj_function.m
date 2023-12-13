function [trajH,trajP,TOL,ALPHA]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,startH,startP,finH,finP,nevol,hostmutationprob)

% This function simulates an evolutionary trajectory from a given starting
% point.

% Set up parameters and vectors to use later:
t_max=100;
tolmin=startH;
tolmax=finH;
alphamin=startP;
alphamax=finP;
res0=101;

% Initial conditions:
strain_totalH = 1;
strain_totalP = 1;
init_pop = [0.1,0.1,0.1,0.1];

% Create vectors of initial conditions:
TolJ = linspace(tolmin,tolmax,res0);
AlphaJ = linspace(alphamin,alphamax,res0);
initialH = find(TolJ>=tolstart,1);
initialP = find(AlphaJ>=alphastart,1);
tol_start = TolJ(initialH);
alpha_start = AlphaJ(initialP);
indexH_start = initialH;
indexP_start = initialP;

% We let the traits evolve for a fixed number of evolutionary timesteps.
% Run the evolutionary simulation:
[~,~,~,~,~,~,~,TOL,ALPHA,~,~] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tolmin,tolmax,tol_start,alphamin,alphamax,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);

% Create output vectors:
trajH=NaN(nevol,1);
trajP=NaN(nevol,1);
for i=1:nevol
    TOLend1=TOL(i,:);
    TOLend=[0 TOLend1 0];
    [~,locsH]=findpeaks(TOLend);
    peaknumH=length(locsH);
    ALPHAend1=ALPHA(i,:);
    ALPHAend=[0 ALPHAend1 0];
    [~,locsP]=findpeaks(ALPHAend);
    peaknumP=length(locsP);
    if peaknumH==1
        trajH(i)=TolJ(locsH-1);
    else
        trajH(i)=NaN;
    end
    if peaknumP==1
        trajP(i)=AlphaJ(locsP-1);
    else
        trajP(i)=NaN;
    end
end


end
