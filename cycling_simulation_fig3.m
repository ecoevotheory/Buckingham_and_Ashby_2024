% This code draws the evolutionary trajectories of juvenile tolerance and
% virulence (and also shows how host and pathogen population densities 
% change with time). 

%% Define Parameters:

rng(10)

% Parameter values:
q=1;
f=1;
bJ=1/5;
g=1;
bA=1/5;
tA=0;
beta0=10;
rJ=0;
rA=0;
a0=5;
c1=0.25;
c2=3;
gamma=0;

tolJmin=0;
tolJmax=1;
tolJstart=0.9;
alphamin=0;
alphamax=10;
alphastart=3.7;

t_max=100;
res0=101;
nevol=10000;

% Set up vectors:
TOLJ=[];
ALPHAJ=[];
DISPREV=[];
N=[];

%% Set up Initial Conditions for a Monomorphic Population

% There is one strain of host and one strain of parasite initially:
strain_totalH = 1;
strain_totalP = 1;

% Initial population distribution (SJ, SA, IJ, IA):
init_pop = [0.1,0.1,0.1,0.1];

% Initial conditions
TolJ = linspace(tolJmin,tolJmax,res0);
AlphaJ = linspace(alphamin,alphamax,res0);
initialH = find(TolJ>=tolJstart,1);
initialP = find(AlphaJ>=alphastart,1);
tolJ_start = TolJ(initialH);
alpha_start = AlphaJ(initialP);
indexH_start = initialH;
indexP_start = initialP;

%% Allow both traits to evolve
% This section can be run multiple times to extend the number of
% evolutionary timesteps without having to start the simulation again. 

% The parasite mutates 20 times faster than the host:
hostmutationprob=1/21;
[tolJ_start,alpha_start,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,TOLJnew,ALPHAJnew,DISPREVnew,Nnew] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tolJmin,tolJmax,tolJ_start,alphamin,alphamax,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
TOLJ=[TOLJ;TOLJnew];
ALPHAJ=[ALPHAJ;ALPHAJnew];
DISPREV=[DISPREV;DISPREVnew];
N=[N;Nnew];

%% Make the Plot

% Plot the host juvenile tolerance trajectory:
clf
aa0=1;
TOLJ0=log10(TOLJ);
TOLJ0(TOLJ0<-aa0)=-aa0;
TOLJ0=(TOLJ0+aa0)/aa0;
subplot(1,5,1)
imagesc(TOLJ0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
ylabel('Evolutionary time','interpreter','latex')
xlabel('Juvenile tolerance, $t_J$','interpreter','latex')
set(gca,'xtick',1:1:res0,'xticklabel',round(TolJ(1:1:res0)*10000)/10000);
set(gca,'xtick',1:(res0-1)/2:res0,'xticklabel',round(TolJ(1:(res0-1)/2:res0)*10000)/10000);
set(gca,'ytick',[0,2000,4000,6000,8000,10000])
title('A')

% Plot the pathogen virulence trajectory:
ab0=1;
ALPHAJ0=log10(ALPHAJ);
ALPHAJ0(ALPHAJ0<-ab0)=-ab0;
ALPHAJ0=(ALPHAJ0+ab0)/ab0;
subplot(1,5,2)
imagesc(ALPHAJ0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
xlabel('Virulence, $\alpha$','interpreter','latex')
set(gca,'xtick',1:1:res0,'xticklabel',round(AlphaJ(1:1:res0)*10000)/10000);
set(gca,'xtick',1:(res0-1)/2:res0,'xticklabel',round(AlphaJ(1:(res0-1)/2:res0)*10000)/10000);
set(gca,'ytick',[0,2000,4000,6000,8000,10000])
title('B')

% Plot the host and pathogen densities:
red=1/255*[215,48,39];
blue=1/255*[69,117,180];
pathogen_density=DISPREV.*N;
subplot(1,5,3)
plot(pathogen_density,linspace(1,length(pathogen_density),length(pathogen_density)),'color',red,'linewidth',2)
hold on
plot(N,linspace(1,length(N),length(N)),'color',blue,'linewidth',2)
xlabel('Population density','interpreter','latex')
xlim([0,0.25])
set(gca,'xtick',[0,0.25])
set(gca,'ytick',[0,2000,4000,6000,8000,10000])
title('C')
