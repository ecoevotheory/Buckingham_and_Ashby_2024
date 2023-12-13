% This code finds and plots the co-singular strategies and their stability 
% for different values of the lifespan (with pathogen baseline 
% transmissibility varying so as to fix disease prevalence). The disease 
% prevalence is also plotted. 

% Define parameter values:
q=1;
f=1;
g=1;
tA=0;
rJ=0;
rA=0;
a0=1;
c1=1;
c2=4;
gamma=1;
maxsingstrats=100;
finP=10;
initvec=[0.1,0.1,0.1,0.1];
orig_tmax=100;

lifespan_min=0.5;
lifespan_max=10;
lifespan_numberofpoints=20;
lifespanvec=linspace(lifespan_min,lifespan_max,lifespan_numberofpoints);

beta0_numberofpoints=10;

% For each value of the lifespan, we want to choose a value of beta0 so
% that the disease prevalence is approximately equal to this target value:
target=0.11;

%% Step 1: Find the co-singular strategies for the lifelong tolerance scenario

LL_tol_vec=NaN(lifespan_numberofpoints,1);
LL_vir_vec=NaN(lifespan_numberofpoints,1);
LL_beta0_vec=NaN(lifespan_numberofpoints,1);
LL_disprev_vec=NaN(lifespan_numberofpoints,1);
% For each value of lifespan, find the co-singular strategies and the value
% of beta0 which makes the disease prevalence approximately equal to the
% target:
for i=1:lifespan_numberofpoints
    disp("Lifespan number " +i)
    bJ=1/lifespanvec(i);
    bA=1/lifespanvec(i);
        
    upper_limit=32;
    lower_limit=0;
    j=1;
    
    % Set up vectors to use later:
    LL_tol_matrix=NaN(1,beta0_numberofpoints);
    LL_vir_matrix=NaN(1,beta0_numberofpoints);
    LL_disprev_matrix=NaN(1,beta0_numberofpoints);
    beta0_matrix=NaN(1,beta0_numberofpoints);
    
    % We "zoom in" on values of beta0 which produce a disease prevalence
    % close to our target:
    while j<beta0_numberofpoints+1
        disp(j)
        beta0=(upper_limit+lower_limit)/2;
        
        % Find co-singular strategies:
        [LL_tol,LL_vir,LL_Hclass,LL_Pclass]=LL_singstrat_function(q,f,bJ,bA,g,beta0,rJ,rA,a0,c1,c2,gamma,finP);
        
        % Where more than one co-singular strategy is found for these
        % parameters, it is always just one co-singular strategy:
        if length(LL_tol)>1
            LL_tol=LL_tol(1);
            LL_vir=LL_vir(1);
            LL_Hclass=LL_Hclass(1);
            LL_Pclass=LL_Pclass(1);
        end
        
        % Find the disease prevalence: 
        if length(LL_tol)==1 && ~isnan(LL_vir) && LL_vir~=inf
            LL_tol_matrix(1,j)=LL_tol;
            LL_vir_matrix(1,j)=LL_vir;
            tJ=LL_tol;
            alpha=LL_vir;
            a=a0*(1-(c1*(1-exp(c2*tJ)))/(1-exp(c2)));
            [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tJ,tJ,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
            LL_disprev=(IJval+IAval)/(SJval+SAval+IJval+IAval);
            LL_disprev_matrix(1,j)=LL_disprev;
            beta0_matrix(1,j)=beta0;
            
            % Adjust the value of beta0 so that the disease prevalence is
            % closer to its target value: 
            if LL_disprev<target
                lower_limit=beta0;
            elseif LL_disprev>target
                upper_limit=beta0;
            end
        end
        
        if isempty(LL_tol) || isnan(LL_vir) || LL_vir==inf
            lower_limit=beta0;
        end
        j=j+1;
        
    end
    
    % Create output vectors:
    if ismember(lower_limit,beta0_matrix) && ismember(upper_limit,beta0_matrix)
        if abs(target-LL_disprev_matrix(beta0_matrix==lower_limit))<abs(LL_disprev_matrix(beta0_matrix==upper_limit)-target)
            LL_tol_vec(i)=LL_tol_matrix(beta0_matrix==lower_limit);
            LL_vir_vec(i)=LL_vir_matrix(beta0_matrix==lower_limit);
            LL_beta0_vec(i)=beta0_matrix(beta0_matrix==lower_limit);
            LL_disprev_vec(i)=LL_disprev_matrix(beta0_matrix==lower_limit);
        else
            LL_tol_vec(i)=LL_tol_matrix(beta0_matrix==upper_limit);
            LL_vir_vec(i)=LL_vir_matrix(beta0_matrix==upper_limit);
            LL_beta0_vec(i)=beta0_matrix(beta0_matrix==upper_limit);
            LL_disprev_vec(i)=LL_disprev_matrix(beta0_matrix==upper_limit);
        end
    end
    
end


%% Step 2: Find the co-singular strategies for the juvenile tolerance scenario

% Set up vectors to use later:
JL_tol_vec=NaN(lifespan_numberofpoints,1);
JL_vir_vec=NaN(lifespan_numberofpoints,1);
JL_beta0_vec=NaN(lifespan_numberofpoints,1);
JL_disprev_vec=NaN(lifespan_numberofpoints,1);

% For each value of lifespan, find the co-singular strategies and the value
% of beta0 which makes the disease prevalence approximately equal to the
% target:
for i=1:lifespan_numberofpoints
    disp("Lifespan number " +i)
    bJ=1/lifespanvec(i);
    bA=1/lifespanvec(i);
        
    upper_limit=32;
    lower_limit=0;
    j=1;
    
    JL_tol_matrix=NaN(1,beta0_numberofpoints);
    JL_vir_matrix=NaN(1,beta0_numberofpoints);
    JL_disprev_matrix=NaN(1,beta0_numberofpoints);
    beta0_matrix=NaN(1,beta0_numberofpoints);
    
    % We "zoom in" on values of beta0 which make the disease prevalence
    % closer to the target:
    while j<beta0_numberofpoints+1
        disp(j)
        beta0=(upper_limit+lower_limit)/2;
        
        % Find the co-singular strategies:
        [JL_tol,JL_vir,JL_Hclass,JL_Pclass]=JL_singstrat_function(q,f,bJ,bA,g,tA,beta0,rJ,rA,a0,c1,c2,gamma,finP);
        
        % Where more than one co-singular strategy is found for these
        % parameters, it is always just one co-singular strategy:
        if length(JL_tol)>1
            JL_tol=JL_tol(1);
            JL_vir=JL_vir(1);
            JL_Hclass=JL_Hclass(1);
            JL_Pclass=JL_Pclass(1);
        end
        
        % Determine the disease prevalence:
        if length(JL_tol)==1 && ~isnan(JL_vir) && JL_vir~=inf
            JL_tol_matrix(1,j)=JL_tol;
            JL_vir_matrix(1,j)=JL_vir;
            tJ=JL_tol;
            alpha=JL_vir;
            a=a0*(1-(c1*(1-exp(c2*tJ)))/(1-exp(c2)));
            [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tJ,tA,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
            JL_disprev=(IJval+IAval)/(SJval+SAval+IJval+IAval);
            JL_disprev_matrix(1,j)=JL_disprev;
            beta0_matrix(1,j)=beta0;
            
            % Adjust the value of beta0 to give a disease prevalence closer
            % to the target:
            if JL_disprev<target
                lower_limit=beta0;
            elseif JL_disprev>target
                upper_limit=beta0;
            end
        end
        
        if isempty(JL_tol) || isnan(JL_vir) || JL_vir==inf
            lower_limit=beta0;
        end
        j=j+1;
        
    end
    
    % Create output vectors:
    if ismember(lower_limit,beta0_matrix) && ismember(upper_limit,beta0_matrix)
        if abs(target-JL_disprev_matrix(beta0_matrix==lower_limit))<abs(JL_disprev_matrix(beta0_matrix==upper_limit)-target)
            JL_tol_vec(i)=JL_tol_matrix(beta0_matrix==lower_limit);
            JL_vir_vec(i)=JL_vir_matrix(beta0_matrix==lower_limit);
            JL_beta0_vec(i)=beta0_matrix(beta0_matrix==lower_limit);
            JL_disprev_vec(i)=JL_disprev_matrix(beta0_matrix==lower_limit);
        else
            JL_tol_vec(i)=JL_tol_matrix(beta0_matrix==upper_limit);
            JL_vir_vec(i)=JL_vir_matrix(beta0_matrix==upper_limit);
            JL_beta0_vec(i)=beta0_matrix(beta0_matrix==upper_limit);
            JL_disprev_vec(i)=JL_disprev_matrix(beta0_matrix==upper_limit);
        end
    end
    
end


%% Section 3 - Create plot

% Colours for plotting:
blue=1/255*[69,117,180];
red=1/255*[215,48,39];

fig=figure;
set(fig,'defaultAxesColorOrder',[0,0,0;0,0,0])

% Plot host tolerance
subplot(2,1,1)
hold on
plot(lifespanvec,LL_tol_vec,'linewidth',3,'linestyle','-','color',blue)
plot(lifespanvec,JL_tol_vec,'linewidth',3,'linestyle','-','color',red)

xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Tolerance, $\tau_i$','interpreter','latex')
xlim([0,max(lifespanvec)])
ylim([0,1])
set(gca,'xtick',[0,1,2,3,4,5,6,7,8,9,10])
text(0.2,0.92,"B (i)",'fontsize',20)
pbaspect([2,1,1])
box on

% Plot where the lifespan gets too low to sustain the host population:
lifespan_min=2/(-g+sqrt(g^2+4*a0*g));
lifespan_min_plotter=lifespanvec(find(lifespanvec>lifespan_min,1));
plot(lifespan_min_plotter+zeros(100,1),linspace(0,1,100),':k','linewidth',2)

% Plot disease prevalence:
yyaxis right
plot(lifespanvec,LL_disprev_vec,'linewidth',2,'linestyle','--','color',blue)
plot(lifespanvec,JL_disprev_vec,'linewidth',2,'linestyle',':','color',red)
ylim([0,1])
ylabel('Disease prevalence')

% Plot parasite virulence
subplot(2,1,2)
hold on
plot(lifespanvec,LL_vir_vec,'linewidth',3,'linestyle','-','color',blue)
plot(lifespanvec,JL_vir_vec,'linewidth',3,'linestyle','-','color',red)

xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Virulence, $\alpha$','interpreter','latex')
ylim([0,5])
xlim([0,max(lifespanvec)])
set(gca,'xtick',[0,1,2,3,4,5,6,7,8,9,10])
text(0.2,4.6,"B (ii)",'fontsize',20)
pbaspect([2,1,1])
box on

% Plot where the lifespan gets too low to sustain the host population:
lifespan_min=2/(-g+sqrt(g^2+4*a0*g));
lifespan_min_plotter=lifespanvec(find(lifespanvec>lifespan_min,1));
plot(lifespan_min_plotter+zeros(100,1),linspace(0,finP,100),':k','linewidth',2)

% Plot disease prevalence:
yyaxis right
plot(lifespanvec,LL_disprev_vec,'linewidth',2,'linestyle','--','color',blue)
plot(lifespanvec,JL_disprev_vec,'linewidth',2,'linestyle',':','color',red)
ylim([0,1])
ylabel('Disease prevalence')
