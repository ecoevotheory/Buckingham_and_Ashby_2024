% This code uses simulations to determine whether cycling occurs for
% different values of the host lifespan and then creates a plot showing the
% maximum and minimum extent of the cycling trait for each value of the
% host lifespan. 

%% Section 1

% First, we will find where cycling occurs, along with the co-CSS's which
% occur via a bifurcation as lifespan decreases (the cycles reduce in
% amplitude until they no longer cycle but form CSS's). 

% Set up parameter values:
lifespanvec=linspace(4,10,25);
q=1;
f=1;
g=1;
tA=0;
beta0=10;
rJ=0;
rA=0;
a0=5;
c1=0.275;
c2=3;
gamma=0;
scenario=2; % Juvenile tolerance model

tolJmin=0;
tolJmax=1;
tolJstart=0.9;
alphamin=0;
alphamax=10;
alphastart=3.7;

t_max=100;
res0=101;
nevol=5000;
hostmutationprob=1/21;
init_pop = [0.1,0.1,0.1,0.1];

% Set up vectors to use later:
cycle_min_vec=NaN(length(lifespanvec),1);
cycle_max_vec=NaN(length(lifespanvec),1);
cycle_both_vec=NaN(length(lifespanvec),1);
tol_min_vec=NaN(length(lifespanvec),1);
tol_max_vec=NaN(length(lifespanvec),1);
tol_both_vec=NaN(length(lifespanvec),1);

% For each value of the lifespan, determine whether or not the simulation
% cycles:
parfor j=1:length(lifespanvec)
    bJ=1/lifespanvec(j);
    bA=1/lifespanvec(j);
    
    disp(j)

    % Set up Initial Conditions
    strain_totalH = 1;
    strain_totalP = 1;
    TolJ = linspace(tolJmin,tolJmax,res0);
    Alpha = linspace(alphamin,alphamax,res0);
    initialH = find(TolJ>=tolJstart,1);
    initialP = find(Alpha>=alphastart,1);
    tolJ_start = TolJ(initialH);
    alpha_start = Alpha(initialP);
    indexH_start = initialH;
    indexP_start = initialP;
    
    % Allow both traits to evolve
    rng(1,'twister')
    [tolJ_start,alpha_start,~,strain_totalH,strain_totalP,indexH_start,indexP_start,TOLJ,ALPHA,DISPREV,N] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tolJmin,tolJmax,tolJ_start,alphamin,alphamax,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
    
    % Find the virulence trait value at different points in time
    traitvec=NaN(nevol,1);
    for i=1:nevol
        ALPHAend1=ALPHA(i,:);
        ALPHAend=[0 ALPHAend1 0];
        [~,SSlocs]=findpeaks(ALPHAend);
        singstrats=Alpha(SSlocs-1);
        % If there is more than one end-point then it must be a branching
        % point. This is denoted by "-1".
        if length(SSlocs)>1
            traitvec(i)=-1;
            % If there is one end-point then this is the current value of the
            % evolving trait:
        elseif length(SSlocs)==1
            traitvec(i)=singstrats;
            % Otherwise, the disease has gone extinct. This is denoted by "-2".
        elseif isempty(SSlocs)
            traitvec(i)=-2;
        end
    end
    levelout=0;
    traitvec_end=traitvec(end-100:end);
    if sum(traitvec_end==traitvec_end(end))==length(traitvec_end)
        levelout=1;
    end
    
    % Determine whether or not the simulation is cycling:
    cycling=0;
    if sum(traitvec==-1)<500 && sum(traitvec==-2)==0 && levelout==0
        truncated_traitvec=traitvec;
        truncated_traitvec(1:1500)=[];
        trait_start=traitvec(1000);
        if (max(truncated_traitvec)-trait_start>0.2 || trait_start-min(truncated_traitvec)>0.2) && sum(truncated_traitvec==trait_start)~=0
            cycling=1;
        end
    else
        truncated_traitvec=traitvec;
        truncated_traitvec(1:1500)=[];
    end
    
    % Now determine whether tolerance branches or reaches a stable
    % equilibrium: 
    tolvec=NaN(nevol,1);
    for i=1:nevol
        TOLJend1=TOLJ(i,:);
        TOLJend=[0 TOLJend1 0];
        [~,SSlocs]=findpeaks(TOLJend);
        singstrats=TolJ(SSlocs-1);
        if length(SSlocs)>1
            tolvec(i)=-1;
            % If there is one end-point then this is the current value of the
            % evolving trait:
        elseif length(SSlocs)==1
            tolvec(i)=singstrats;
            % Otherwise, the host has gone extinct. This is denoted by "-2".
        elseif isempty(SSlocs)
            tolvec(i)=-2;
        end
    end
    truncated_tolvec=tolvec;
    truncated_tolvec(1:1500)=[];
    
    % If there is cycling, determine the maximum and minimum extents of the
    % cyles:
    if cycling==1
        truncated_traitvec(truncated_traitvec==-1)=[];
        cycle_min_vec(j)=min(truncated_traitvec);
        cycle_max_vec(j)=max(truncated_traitvec);
        
        truncated_tolvec(truncated_tolvec==-1)=[];
        tol_min_vec(j)=min(truncated_tolvec);
        tol_max_vec(j)=max(truncated_tolvec);
    elseif traitvec(end)==10
        cycle_both_vec(j)=10;
        tol_both_vec(j)=1;
    else
        cycle_min_vec(j)=traitvec(end);
        cycle_max_vec(j)=traitvec(end);
        cycle_both_vec(j)=traitvec(end);
        
        tol_min_vec(j)=tolvec(end);
        tol_max_vec(j)=tolvec(end);
        tol_both_vec(j)=tolvec(end);
    end

end

% Once the curves reach their maximum values, they stay there:
for j=1:length(lifespanvec)-1
    if cycle_both_vec(j)==10
        for k=j+1:length(lifespanvec)
            cycle_both_vec(k)=10;
            cycle_min_vec(k)=NaN;
            cycle_max_vec(k)=NaN;
            
            tol_both_vec(k)=1;
            tol_min_vec(k)=NaN;
            tol_max_vec(k)=NaN;
        end
    end
end

% Smooth parts of the curves:
smooth_traj1=smoothdata(cycle_max_vec,'movmean',5);
smooth_traj2=smoothdata(cycle_min_vec,'movmean',5);
smooth_traj1=[cycle_max_vec(1:4);smooth_traj1(5:16);cycle_max_vec(17:end)];
smooth_traj2=[cycle_min_vec(1:4);smooth_traj2(5:16);cycle_min_vec(17:end)];
% This point is added to the end of the vector to make the curves join up:
cycle_both_vec(5)=3.6;

% Create final vectors for plotting:
vir_max_vec=smooth_traj1;
vir_min_vec=smooth_traj2;
vir_both_vec=cycle_both_vec;

% Smooth parts of the curves:
smooth_traj3=smoothdata(tol_max_vec,'movmean',5);
smooth_traj4=smoothdata(tol_min_vec,'movmean',5);
smooth_traj3=[tol_max_vec(1:4);smooth_traj3(5:16);tol_max_vec(17:end)];
smooth_traj4=[tol_min_vec(1:4);smooth_traj4(5:16);tol_min_vec(17:end)];
% This point is added to the end of the vector to make the curves join up:
tol_both_vec(5)=0.895;

% Create final vectors for plotting:
tol_max_vec=smooth_traj3;
tol_min_vec=smooth_traj4;

%% Section 2

% Now we will find the co-CSS's which lie below the cycles.

tol_min=tolJmin;
tol_max=tolJmax;
alpha_min=alphamin;
finished=0;
eps=0.02;
maxsingstrats=100;

tol_matrix=NaN(25,100);
vir_matrix=NaN(25,100);
Hclass_matrix=NaN(25,100);
Pclass_matrix=NaN(25,100);

for steps=1:25

    disp("steps = "+steps)
    
    bJ=1/lifespanvec(steps);
    bA=1/lifespanvec(steps);
    
    if gamma==0
        finP=10;
    elseif gamma==1
        finP=25;
    end
    alpha_max=finP;
    
    % Section 1 - Find co-singular strategies and their stability
    
    % The classifications represent:
    % 1: a co-CSS (evo stab, conv stab, strong conv stab)
    % 2: a repeller (conv unstab, strong conv unstab)
    % 3: a branching point (evo unstab, conv stab, strong conv stab)
    % 4: an evolutionarily stable ambiguous case (evo stab, conv stab
    % when parasite evolves much faster than host, strong conv unstab)
    % 5: an evolutionarily unstable ambiguous case (evo unstab, conv
    % stab when parasite evolves much faster than host, strong conv unstab)
    % 10n+1: the classification "n" reached from analytical calculations with
    % simulations revealing cycling
    % 10n+2: the classification "n" reached from analytical calculations with
    % simulations revealing branching
    % 10n+3: the classification "n" reached from analytical calculations with
    % simulations revealing an attractor
    % 10n+4: the classification "n" reached from analytical calculations with
    % simulations revealing a repeller
    % 10n+5: the classification "n" reached from analytical calculations with
    % simulations revealing parasite extinction
    
    % Find the co-singular strategies for tolerance and virulence:
    if scenario==1
        [tol,vir,Hclass,Pclass]=LL_singstrat_function(q,f,bJ,bA,g,beta0,rJ,rA,a0,c1,c2,gamma,finP);
    elseif scenario==2
        [tol,vir,Hclass,Pclass]=JL_singstrat_function(q,f,bJ,bA,g,tA,beta0,rJ,rA,a0,c1,c2,gamma,finP);
    end
    
    % Section 2 - Determine whether cycling occurs around any of the co-singular strategies
    
    % Set up initial conditions:
    strain_totalH = 1;
    strain_totalP = 1;
    Tol = linspace(tol_min,tol_max,res0);
    Alpha = linspace(alpha_min,alpha_max,res0);
    
    parfor j=1:length(tol)
        
        if vir(j)~=inf
            
            tolstart=min(tol(j)+eps,tol_max);
            alphastart=min(vir(j)+eps,alpha_max);
            
            % Set up Initial Conditions
            initialH = find(Tol>=tolstart,1);
            initialP = find(Alpha>=alphastart,1);
            tol_start = Tol(initialH);
            alpha_start = Alpha(initialP);
            indexH_start = initialH;
            indexP_start = initialP;
            
            % Allow both traits to evolve
            rng(2)
            ALPHA=[];
            if scenario==1
                [~,~,~,~,~,~,~,TOL,ALPHA,~,~] = LL_simulation_function(t_max,a0,g,q,beta0,c1,c2,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
            elseif scenario==2
                [~,~,~,~,~,~,~,TOL,ALPHA,~,~] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
            end
            
            % Find the virulence trait value at different points in time
            virvec=NaN(nevol,1);
            for i=1:nevol
                ALPHAend1=ALPHA(i,:);
                ALPHAend=[0 ALPHAend1 0];
                [~,SSlocs]=findpeaks(ALPHAend);
                singstrats=Alpha(SSlocs-1);
                % If there is more than one end-point then it must be a branching
                % point. This is denoted by "-1".
                if length(SSlocs)>1
                    virvec(i)=-1;
                    % If there is one end-point then this is the current value of the
                    % evolving trait:
                elseif length(SSlocs)==1
                    virvec(i)=singstrats;
                    % Otherwise, the disease has gone extinct. This is denoted by "-2".
                elseif isempty(SSlocs)
                    virvec(i)=-2;
                end
            end
            levelout=0;
            virvec_end=virvec(end-100:end);
            if sum(virvec_end==virvec_end(end))==length(virvec_end)
                levelout=1;
            end
            
            % Determine whether or not virulence is cycling:
            vir_cycling=0;
            if sum(virvec==-1)<500 && sum(virvec==-2)==0 && levelout==0
                truncated_virvec=virvec;
                truncated_virvec(1:1500)=[];
                vir_start=virvec(1000);
                if (max(truncated_virvec)-vir_start>0.2 || vir_start-min(truncated_virvec)>0.2) && sum(truncated_virvec==vir_start)~=0
                    vir_cycling=1;
                end
            end
            
            % Find the tolerance trait value at different points in time
            tolvec=NaN(nevol,1);
            for i=1:nevol
                TOLend1=TOL(i,:);
                TOLend=[0 TOLend1 0];
                [~,SSlocs]=findpeaks(TOLend);
                singstrats=Tol(SSlocs-1);
                if length(SSlocs)>1
                    tolvec(i)=-1;
                    % If there is one end-point then this is the current value of the
                    % evolving trait:
                elseif length(SSlocs)==1
                    tolvec(i)=singstrats;
                    % Otherwise, the host has gone extinct. This is denoted by "-2".
                elseif isempty(SSlocs)
                    tolvec(i)=-2;
                end
            end
            
            % Determine whether or not tolerance is cycling:
            tol_cycling=0;
            if sum(tolvec==-1)<500 && sum(tolvec==-2)==0 && levelout==0
                truncated_tolvec=tolvec;
                truncated_tolvec(1:1500)=[];
                tol_start=tolvec(1000);
                if (max(truncated_tolvec)-tol_start>0.2 || tol_start-min(truncated_tolvec)>0.2) && sum(truncated_tolvec==tol_start)~=0
                    tol_cycling=1;
                end
            end
            
            if vir_cycling==1
                Hclass(j)=Hclass(j)*10+1;
            end
            if tol_cycling==1
                Pclass(j)=Pclass(j)*10+1;
            end
            
        end
    end
    
    % Section 3 - Investigate ambiguous convergence stability cases
    
    parfor j=1:length(tol)
        if (Hclass(j)==4 || Pclass(j)==4 || Hclass(j)==5 || Pclass(j)==5) && Hclass(j)~=21 && Pclass(j)~=21 && Hclass(j)~=41 && Pclass(j)~=41 && Hclass(j)~=51 && Pclass(j)~=51
            
            tolstart=min(tol(j)+((tol_max-tol_min)/res0)*10,tol_max);
            alphastart=min(vir(j)+((alpha_max-alpha_min)/res0)*10,alpha_max);
            
            % Set up Initial Conditions
            initialH = find(Tol>=tolstart,1);
            initialP = find(Alpha>=alphastart,1);
            tol_start = Tol(initialH);
            alpha_start = Alpha(initialP);
            indexH_start = initialH;
            indexP_start = initialP;
            
            % Allow both traits to evolve
            rng(1)
            if scenario==1
                [~,~,~,~,~,~,~,TOL1,ALPHA1,~,~] = LL_simulation_function(t_max,a0,g,q,beta0,c1,c2,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
            elseif scenario==2
                [~,~,~,~,~,~,~,TOL1,ALPHA1,~,~] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
            end
            
            % Find whether the virulence trait has evolved towards or away from
            % the co-singular strategy:
            ALPHAend1=ALPHA1(nevol,:);
            ALPHAend=[0 ALPHAend1 0];
            [~,SSlocs]=findpeaks(ALPHAend);
            if length(SSlocs)>1 % suggests branching
                Pclass(j)=10*Pclass(j)+2;
            elseif length(SSlocs)==1
                if abs(SSlocs-vir(j))<abs(alphastart-vir(j))
                    Pclass(j)=10*Pclass(j)+3;
                else
                    Pclass(j)=10*Pclass(j)+4;
                end
            elseif isempty(SSlocs) % suggests parasite extinction
                Pclass(j)=10*Pclass(j)+5;
            end
            
            % Find whether the tolerance trait has evolved towards or away from
            % the co-singular strategy:
            TOLend1=TOL1(nevol,:);
            TOLend=[0 TOLend1 0];
            [~,SSlocs]=findpeaks(TOLend);
            if length(SSlocs)>1 % suggests branching
                Hclass(j)=10*Hclass(j)+2;
            elseif length(SSlocs)==1
                if abs(SSlocs-tol(j))<abs(tolstart-tol(j))
                    Hclass(j)=10*Hclass(j)+3;
                else
                    Hclass(j)=10*Hclass(j)+4;
                end
            elseif isempty(SSlocs) % suggests parasite extinction
                Hclass(j)=10*Hclass(j)+5;
            end
            
        end
    end
    
    % Section 4 - Output final classification
    
    final_class=NaN(length(tol),1);
    for j=1:length(tol)
        if sum(vir==inf)==length(vir)
            final_class(j)=[]; % extinction
        elseif Hclass(j)==Pclass(j) && ismember(Hclass(j),[1,43])
            final_class(j)=1; % co-CSS
        elseif Hclass(j)==Pclass(j) && ismember(Hclass(j),[2,44,54])
            final_class(j)=2; % repeller
        elseif ismember(Hclass(j),[3,42,52,53]) || ismember(Pclass(j),[3,42,52,53])
            final_class(j)=3; % branching point
        elseif ismember(Hclass(j),[11,21,31,41,51]) || ismember(Pclass(j),[11,21,31,41,51])
            final_class(j)=4; % cycling
        elseif Hclass(j)==Pclass(j) && ismember(Hclass(j),[45,55])
            final_class(j)=[]; % host extinction
        elseif ismember(Pclass(j),[45,55])
            final_class(j)=[]; % parasite extinction
        elseif ismember(Hclass(j),[2,24,44,54]) || ismember(Pclass(j),[2,24,44,54])
            final_class(j)=2; % repeller
        else
            disp("possible error")
            final_class(j)=[];
        end
    end
    
    if isempty(final_class)
        classification=1; % EXTINCTION
    end
    if ismember(3,final_class) && ismember(4,final_class)
        classification=234;
    elseif ismember(2,final_class) && ismember(4,final_class)
        classification=24;
    elseif ismember(2,final_class) && ismember(3,final_class)
        classification=34;
    elseif ismember(4,final_class) && ismember(1,final_class)
        classification=24;
    elseif ismember(3,final_class) && ismember(1,final_class)
        classification=34;
    elseif ismember(2,final_class)
        classification=4;
    elseif ismember(4,final_class)
        classification=2;
    elseif ismember(3,final_class)
        classification=3;
    elseif ismember(1,final_class)
        classification=5;
    end

    tol_matrix(steps,1:length(tol))=tol;
    vir_matrix(steps,1:length(vir))=vir;
    Hclass_matrix(steps,1:length(Hclass))=transpose(Hclass);
    Pclass_matrix(steps,1:length(Pclass))=transpose(Pclass);

end

tol_lowCSS_vec=tol(find(Pclass==1,1));
vir_lowCSS_vec=vir(find(Pclass==1,1));

%% Section 3

% Plot the figure.

% Colours for plotting:
red=1/255*[215,48,39];
orange=1/255*[253,174,97];
blue=1/255*[69,117,180];

% Smooth parts of the curves:
smooth_traj5=smoothdata(vir_lowCSS_vec,'movmean',5);

% Make the plot for virulence:
subplot(2,1,2)
plot(lifespanvec,vir_max_vec,'color',red,'linewidth',3)
hold on
plot(lifespanvec,vir_min_vec,'color',orange,'linewidth',3)
hold on
plot(lifespanvec,vir_both_vec,'color',blue,'linewidth',3)
hold on
plot(lifespanvec,smooth_traj5,'color',blue,'linewidth',3)
xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Virulence, $\alpha$','interpreter','latex')
ylim([0,10])
xlim([4,10])
text(9.6,9,"B",'fontsize',24)
set(gca,'xtick',[4,5,6,7,8,9,10],'fontsize',14)
set(gca,'ytick',[0,5,10])
pbaspect([2,1,1])

% Smooth parts of the curves:
smooth_traj6=smoothdata(tol_lowCSS_vec,'movmean',5);

% Make the plot for tolerance:
subplot(2,1,1)
plot(lifespanvec,tol_max_vec,'color',red,'linewidth',3)
hold on
plot(lifespanvec,tol_min_vec,'color',orange,'linewidth',3)
hold on
plot(lifespanvec,tol_both_vec,'color',blue,'linewidth',3)
hold on
plot(lifespanvec,smooth_traj6,'color',blue,'linewidth',3)
xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Juvenile tolerance, $\tau_J$','interpreter','latex')
ylim([0,1])
xlim([4,10])
text(9.6,0.9,"A",'fontsize',24)
set(gca,'ytick',[0,0.5,1])
set(gca,'xtick',[4,5,6,7,8,9,10],'fontsize',14)
pbaspect([2,1,1])
