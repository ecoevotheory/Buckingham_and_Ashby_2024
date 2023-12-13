% This code finds and plots the co-singular strategies and their stability 
% for different values of the lifespan (with pathogen baseline 
% transmissibility fixed). It also finds and plots the disease prevalence. 

%% Section 1 - Find co-singular strategies and their stability

% Define parameter values:
param_min=0.5;
param_max=10;
numberofpoints=20;
paramvec=linspace(param_min,param_max,numberofpoints);
numberofpoints=length(paramvec);

q=1;
f=1;
g=1;
bJ=1;
bA=1;
beta0=10;
tA=0;
rJ=0;
rA=0;
a0=1;
c1=1;
c2=4;
gamma=1;

maxsingstrats=100;
finP=5;

% Set up vectors to use later:
LL_tol_vec=NaN(numberofpoints,maxsingstrats);
LL_vir_vec=NaN(numberofpoints,maxsingstrats);
LL_Hclass_vec=NaN(numberofpoints,maxsingstrats);
LL_Pclass_vec=NaN(numberofpoints,maxsingstrats);
JL_tol_vec=NaN(numberofpoints,maxsingstrats);
JL_vir_vec=NaN(numberofpoints,maxsingstrats);
JL_Hclass_vec=NaN(numberofpoints,maxsingstrats);
JL_Pclass_vec=NaN(numberofpoints,maxsingstrats);

% For each value of lifespan, find the co-singular strategies:
parfor i=1:numberofpoints
    disp(i)
    bJ=1/paramvec(i);
    bA=1/paramvec(i);
    
    % Set up vectors to use later:
    LL_tol_vec_temp=NaN(1,maxsingstrats);
    LL_vir_vec_temp=NaN(1,maxsingstrats);
    LL_Hclass_vec_temp=NaN(1,maxsingstrats);
    LL_Pclass_vec_temp=NaN(1,maxsingstrats);
    JL_tol_vec_temp=NaN(1,maxsingstrats);
    JL_vir_vec_temp=NaN(1,maxsingstrats);
    JL_Hclass_vec_temp=NaN(1,maxsingstrats);
    JL_Pclass_vec_temp=NaN(1,maxsingstrats);
    
    % Find the co-singular strategies for lifelong tolerance and virulence:
    [LL_tol,LL_vir,LL_Hclass,LL_Pclass]=LL_singstrat_function(q,f,bJ,bA,g,beta0,rJ,rA,a0,c1,c2,gamma,finP);
    if length(LL_tol)>maxsingstrats
        disp("Too many singular strategies")
    end
    LL_tol_vec_temp(1:length(LL_tol))=LL_tol;
    LL_vir_vec_temp(1:length(LL_vir))=LL_vir;
    LL_Hclass_vec_temp(1:length(LL_Hclass))=LL_Hclass;
    LL_Pclass_vec_temp(1:length(LL_Pclass))=LL_Pclass;
    
    % Find the co-singular strategies for juvenile tolerance and virulence:
    [JL_tol,JL_vir,JL_Hclass,JL_Pclass]=JL_singstrat_function(q,f,bJ,bA,g,tA,beta0,rJ,rA,a0,c1,c2,gamma,finP);
    if length(JL_tol)>maxsingstrats
        disp("Too many singular strategies")
    end
    JL_tol_vec_temp(1:length(JL_tol))=JL_tol;
    JL_vir_vec_temp(1:length(JL_vir))=JL_vir;
    JL_Hclass_vec_temp(1:length(JL_Hclass))=JL_Hclass;
    JL_Pclass_vec_temp(1:length(JL_Pclass))=JL_Pclass;
    
    % Create output vectors:
    LL_tol_vec(i,:)=LL_tol_vec_temp;
    LL_vir_vec(i,:)=LL_vir_vec_temp;
    LL_Hclass_vec(i,:)=LL_Hclass_vec_temp;
    LL_Pclass_vec(i,:)=LL_Pclass_vec_temp;
    JL_tol_vec(i,:)=JL_tol_vec_temp;
    JL_vir_vec(i,:)=JL_vir_vec_temp;
    JL_Hclass_vec(i,:)=JL_Hclass_vec_temp;
    JL_Pclass_vec(i,:)=JL_Pclass_vec_temp;
    
end

%% Section 2 - Format data for plotting

% Set up vectors to use later:
LL_vir_vec(LL_vir_vec==inf)=finP+1;
JL_vir_vec(JL_vir_vec==inf)=finP+1;

LL_tol_CSS_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_tol_rep_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_tol_BP_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_tol_evostabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_tol_evounstabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_vir_CSS_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_vir_rep_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_vir_BP_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_vir_evostabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);
LL_vir_evounstabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);

JL_tol_CSS_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_tol_rep_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_tol_BP_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_tol_evostabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_tol_evounstabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_vir_CSS_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_vir_rep_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_vir_BP_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_vir_evostabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);
JL_vir_evounstabconvamb_plotter=NaN(numberofpoints*maxsingstrats,2);

% Include each co-singular strategy in a vector of CSS's, repellers,
% branching points or singular strategies with ambiguous convergence
% stability:
for i=1:numberofpoints
    for j=1:maxsingstrats
        
        % Lifelong tolerance / virulence - host tolerance
        if LL_Hclass_vec(i,j)==1
            LL_tol_CSS_plotter(i+(j-1)*numberofpoints,1)=LL_tol_vec(i,j);
            LL_tol_CSS_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Hclass_vec(i,j)==2
            LL_tol_rep_plotter(i+(j-1)*numberofpoints,1)=LL_tol_vec(i,j);
            LL_tol_rep_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Hclass_vec(i,j)==3
            LL_tol_BP_plotter(i+(j-1)*numberofpoints,1)=LL_tol_vec(i,j);
            LL_tol_BP_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Hclass_vec(i,j)==4
            LL_tol_evostabconvamb_plotter(i+(j-1)*numberofpoints,1)=LL_tol_vec(i,j);
            LL_tol_evostabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Hclass_vec(i,j)==5
            LL_tol_evounstabconvamb_plotter(i+(j-1)*numberofpoints,1)=LL_tol_vec(i,j);
            LL_tol_evounstabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        end
        
        % Lifelong tolerance / virulence - parasite virulence
        if LL_Pclass_vec(i,j)==1 && LL_vir_vec(i,j)~=finP+1
            LL_vir_CSS_plotter(i+(j-1)*numberofpoints,1)=LL_vir_vec(i,j);
            LL_vir_CSS_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Pclass_vec(i,j)==2
            LL_vir_rep_plotter(i+(j-1)*numberofpoints,1)=LL_vir_vec(i,j);
            LL_vir_rep_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Pclass_vec(i,j)==3
            LL_vir_BP_plotter(i+(j-1)*numberofpoints,1)=LL_vir_vec(i,j);
            LL_vir_BP_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Pclass_vec(i,j)==4
            LL_vir_evostabconvamb_plotter(i+(j-1)*numberofpoints,1)=LL_vir_vec(i,j);
            LL_vir_evostabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_Pclass_vec(i,j)==5
            LL_vir_evounstabconvamb_plotter(i+(j-1)*numberofpoints,1)=LL_vir_vec(i,j);
            LL_vir_evounstabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif LL_vir_vec(i,j)==finP+1
            LL_vir_CSS_plotter(i+(j-1)*numberofpoints,1)=LL_vir_vec(i,j);
            LL_vir_CSS_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
            LL_vir_rep_plotter(i+(j-1)*numberofpoints,1)=0;
            LL_vir_rep_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        end
        
        
        % Juvenile tolerance / virulence - host tolerance
        if JL_Hclass_vec(i,j)==1
            JL_tol_CSS_plotter(i+(j-1)*numberofpoints,1)=JL_tol_vec(i,j);
            JL_tol_CSS_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Hclass_vec(i,j)==2
            JL_tol_rep_plotter(i+(j-1)*numberofpoints,1)=JL_tol_vec(i,j);
            JL_tol_rep_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Hclass_vec(i,j)==3
            JL_tol_BP_plotter(i+(j-1)*numberofpoints,1)=JL_tol_vec(i,j);
            JL_tol_BP_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Hclass_vec(i,j)==4
            JL_tol_evostabconvamb_plotter(i+(j-1)*numberofpoints,1)=JL_tol_vec(i,j);
            JL_tol_evostabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Hclass_vec(i,j)==5
            JL_tol_evounstabconvamb_plotter(i+(j-1)*numberofpoints,1)=JL_tol_vec(i,j);
            JL_tol_evounstabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        end
        
        % Juvenile tolerance / virulence - parasite virulence
        if JL_Pclass_vec(i,j)==1
            JL_vir_CSS_plotter(i+(j-1)*numberofpoints,1)=JL_vir_vec(i,j);
            JL_vir_CSS_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Pclass_vec(i,j)==2
            JL_vir_rep_plotter(i+(j-1)*numberofpoints,1)=JL_vir_vec(i,j);
            JL_vir_rep_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Pclass_vec(i,j)==3
            JL_vir_BP_plotter(i+(j-1)*numberofpoints,1)=JL_vir_vec(i,j);
            JL_vir_BP_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Pclass_vec(i,j)==4
            JL_vir_evostabconvamb_plotter(i+(j-1)*numberofpoints,1)=JL_vir_vec(i,j);
            JL_vir_evostabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        elseif JL_Pclass_vec(i,j)==5
            JL_vir_evounstabconvamb_plotter(i+(j-1)*numberofpoints,1)=JL_vir_vec(i,j);
            JL_vir_evounstabconvamb_plotter(i+(j-1)*numberofpoints,2)=paramvec(i);
        end
        
    end
end

% Where the code has found several co-singular strategies close together,
% take the average of these values:
LL_tol_CSS_plotter2=NaN(numberofpoints,2);
for i=1:numberofpoints
    LL_tol_CSS_plotter2(i,1)=mean(LL_tol_CSS_plotter((LL_tol_CSS_plotter(:,2)==paramvec(i)),1));
    LL_tol_CSS_plotter2(i,2)=paramvec(i);
end
LL_vir_CSS_plotter2=NaN(numberofpoints,2);
for i=1:numberofpoints
    LL_vir_CSS_plotter2(i,1)=mean(LL_vir_CSS_plotter((LL_vir_CSS_plotter(:,2)==paramvec(i)),1));
    LL_vir_CSS_plotter2(i,2)=paramvec(i);
end
LL_tol_CSS_plotter=LL_tol_CSS_plotter2;
LL_vir_CSS_plotter=LL_vir_CSS_plotter2;

%% Section 3 - Find disease prevalence

JL_disprev=NaN(length(JL_tol_CSS_plotter),2);
for i=1:size(JL_tol_CSS_plotter,1)
    
    bJ=1/JL_tol_CSS_plotter(i,2);
    bA=bJ;
    tJ=JL_tol_CSS_plotter(i,1);
    alpha=JL_vir_CSS_plotter(i,1);
    
    initvec=[0.1,0.1,0.1,0.1];
    orig_tmax=100;
    a=a0*(1-(c1*(1-exp(c2*tJ)))/(1-exp(c2)));
    if ~isnan(tJ)
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tJ,tA,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
        
        JL_disprev(i,1)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        JL_disprev(i,2)=1/bJ;
    end
    
end

LL_disprev=NaN(length(LL_tol_CSS_plotter),2);
for i=1:size(LL_tol_CSS_plotter,1)
    
    bJ=1/LL_tol_CSS_plotter(i,2);
    bA=bJ;
    tJ=LL_tol_CSS_plotter(i,1);
    alpha=LL_vir_CSS_plotter(i,1);
    
    initvec=[0.1,0.1,0.1,0.1];
    orig_tmax=100;
    a=a0*(1-(c1*(1-exp(c2*tJ)))/(1-exp(c2)));
    if ~isnan(tJ)
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tJ,tJ,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
        
        LL_disprev(i,1)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        LL_disprev(i,2)=1/bJ;
    end
    
end

%% Section 4 - Create plot

% Colours for plotting:
blue=1/255*[69,117,180];
red=1/255*[215,48,39];

fig=figure;
set(fig,'defaultAxesColorOrder',[0,0,0;0,0,0])

% Plot host tolerance
subplot(2,1,1)
hold on
plot(LL_tol_CSS_plotter(:,2),LL_tol_CSS_plotter(:,1),'linewidth',3,'linestyle','-','color',blue)
plot(LL_tol_rep_plotter(:,2),LL_tol_rep_plotter(:,1),'linewidth',3,'linestyle',':','color',blue)
plot(LL_tol_BP_plotter(:,2),LL_tol_BP_plotter(:,1),'linewidth',3,'linestyle','--','color',blue)
plot(LL_tol_evostabconvamb_plotter(:,2),LL_tol_evostabconvamb_plotter(:,1),'linewidth',1,'linestyle','-','color',blue)
plot(LL_tol_evounstabconvamb_plotter(:,2),LL_tol_evounstabconvamb_plotter(:,1),'linewidth',1,'linestyle',':','color',blue)

plot(JL_tol_CSS_plotter(:,2),JL_tol_CSS_plotter(:,1),'linewidth',3,'linestyle','-','color',red)
plot(JL_tol_rep_plotter(:,2),JL_tol_rep_plotter(:,1),'linewidth',3,'linestyle',':','color',red)
plot(JL_tol_BP_plotter(:,2),JL_tol_BP_plotter(:,1),'linewidth',3,'linestyle','--','color',red)
plot(JL_tol_evostabconvamb_plotter(:,2),JL_tol_evostabconvamb_plotter(:,1),'linewidth',1,'linestyle','-','color',red)
plot(JL_tol_evounstabconvamb_plotter(:,2),JL_tol_evounstabconvamb_plotter(:,1),'linewidth',1,'linestyle',':','color',red)

xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Tolerance, $\tau_i$','interpreter','latex')
xlim([0,max(paramvec)])
ylim([0,1])
set(gca,'xtick',[0,1,2,3,4,5,6,7,8,9,10])
text(0.2,0.92,"A (i)",'fontsize',20)
pbaspect([2,1,1])
box on

% Plot where the lifespan gets too low to sustain the host population:
lifespan_min=2/(-g+sqrt(g^2+4*a0*g));
lifespan_min_plotter=paramvec(find(paramvec>lifespan_min,1));
plot(lifespan_min_plotter+zeros(100,1),linspace(0,1,100),':k','linewidth',2)

yyaxis right
plot(LL_disprev(:,2),LL_disprev(:,1),'linewidth',2,'linestyle','--','color',blue)
plot(JL_disprev(:,2),JL_disprev(:,1),'linewidth',2,'linestyle',':','color',red)
ylim([0,1])
ylabel('Disease prevalence')

% Plot parasite virulence
subplot(2,1,2)
hold on
h1=plot(LL_vir_CSS_plotter(:,2),LL_vir_CSS_plotter(:,1),'linewidth',3,'linestyle','-','color',blue);
plot(LL_vir_rep_plotter(:,2),LL_vir_rep_plotter(:,1),'linewidth',3,'linestyle',':','color',blue)
plot(LL_vir_BP_plotter(:,2),LL_vir_BP_plotter(:,1),'linewidth',3,'linestyle','--','color',blue)
plot(LL_vir_evostabconvamb_plotter(:,2),LL_vir_evostabconvamb_plotter(:,1),'linewidth',1,'linestyle','-','color',blue)
plot(LL_vir_evounstabconvamb_plotter(:,2),LL_vir_evounstabconvamb_plotter(:,1),'linewidth',1,'linestyle',':','color',blue)

h2=plot(JL_vir_CSS_plotter(:,2),JL_vir_CSS_plotter(:,1),'linewidth',3,'linestyle','-','color',red);
plot(JL_vir_rep_plotter(:,2),JL_vir_rep_plotter(:,1),'linewidth',3,'linestyle',':','color',red)
plot(JL_vir_BP_plotter(:,2),JL_vir_BP_plotter(:,1),'linewidth',3,'linestyle','--','color',red)
plot(JL_vir_evostabconvamb_plotter(:,2),JL_vir_evostabconvamb_plotter(:,1),'linewidth',1,'linestyle','-','color',red)
plot(JL_vir_evounstabconvamb_plotter(:,2),JL_vir_evounstabconvamb_plotter(:,1),'linewidth',1,'linestyle',':','color',red)

xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Virulence, $\alpha$','interpreter','latex')
ylim([0,finP])
xlim([0,max(paramvec)])
set(gca,'xtick',[0,1,2,3,4,5,6,7,8,9,10])
text(0.2,4.6,"A (ii)",'fontsize',20)
pbaspect([2,1,1])
box on

% Plot where the lifespan gets too low to sustain the host population:
lifespan_min=2/(-g+sqrt(g^2+4*a0*g));
lifespan_min_plotter=paramvec(find(paramvec>lifespan_min,1));
plot(lifespan_min_plotter+zeros(100,1),linspace(0,finP,100),':k','linewidth',2)

yyaxis right
plot(LL_disprev(:,2),LL_disprev(:,1),'linewidth',2,'linestyle','--','color',blue)
plot(JL_disprev(:,2),JL_disprev(:,1),'linewidth',2,'linestyle',':','color',red)
ylim([0,1])
ylabel('Disease prevalence')
