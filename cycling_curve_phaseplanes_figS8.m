% This code draws a bifurcation diagram showing the effect of lifespan on
% cycling. It also creates three phase planes for different values of
% lifespan. 

%% Set up parameters to use later

% Parameter values:
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
hostmutationprob=1/21;

tolJmin=0;
tolJmax=1;
alphamin=0;
alphamax=10;
t_max=100;

% Colours for plotting:
red=1/255*[215,48,39];
blue=1/255*[69,117,180];
orange=1/255*[217,95,2];
black=1/255*[0,0,0];
grey=1/255*[166,166,166];

%% Bifurcation plots
% For different values of lifespan, determine whether or not the
% co-evolutionary trajectory cycles. If it does, determine the maximum and
% minimum extent of those cycles. Plot these as a bifurcation diagram. 

% Set up parameter values:
lifespanvec=linspace(4,10,25);
tolJstart=0.9;
alphastart=3.7;
res0=101;
nevol=5000;
init_pop = [0.1,0.1,0.1,0.1];

% Set up vectors to use later:
cycle_min_vec=NaN(length(lifespanvec),1);
cycle_max_vec=NaN(length(lifespanvec),1);
cycle_both_vec=NaN(length(lifespanvec),1);
tol_min_vec=NaN(length(lifespanvec),1);
tol_max_vec=NaN(length(lifespanvec),1);
tol_both_vec=NaN(length(lifespanvec),1);

% For each value of lifespan, determine whether or not the evolutionary
% trajectory cycles:
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
    
    % Determine whether or not there is cycling:
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
    
    % Now determine whether tolerance branches or goes to a stable
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
    
    % If there is cycling, find the maximum and minimum extents of the
    % cycles:
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

% Plot the figure.

% Smooth parts of the curves:
smooth_traj5=smoothdata(vir_lowCSS_vec,'movmean',5);

% Make the plot for virulence:
subplot(4,3,[4;6])
plot(lifespanvec,vir_max_vec,'color',red,'linewidth',2)
hold on
plot(lifespanvec,vir_min_vec,'color',orange,'linewidth',2)
hold on
plot(lifespanvec,vir_both_vec,'color',blue,'linewidth',2)
hold on
plot(lifespanvec,smooth_traj5,'color',blue,'linewidth',2)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',14)
ylim([0,10])
xlim([4,9])
text(8.85,8.6,"B",'fontsize',18)
set(gca,'xtick',[4,5,6,7,8,9],'fontsize',14)
xlabel('Lifespan, $1/b$','interpreter','latex','fontsize',14)

% Smooth parts of the curves:
smooth_traj6=smoothdata(tol_lowCSS_vec,'movmean',5);

% Make the plot for tolerance:
subplot(4,3,[1;3])
plot(lifespanvec,tol_max_vec,'color',red,'linewidth',2)
hold on
plot(lifespanvec,tol_min_vec,'color',orange,'linewidth',2)
hold on
plot(lifespanvec,tol_both_vec,'color',blue,'linewidth',2)
hold on
plot(lifespanvec,smooth_traj6,'color',blue,'linewidth',2)
xlabel('Lifespan, $1/b$','interpreter','latex','fontsize',14)
ylabel('Tolerance, $t_J$','interpreter','latex','fontsize',14)
ylim([0,1])
xlim([4,9])
text(8.85,0.96,"A",'fontsize',18)
set(gca,'ytick',[0,0.5,1])
set(gca,'xtick',[4,5,6,7,8,9],'fontsize',14)

%% Set up parameters to use for phase planes

initvec=[0.1,0.1,0.01,0.01];
SSres=500;

tolstart=0.95;
alphastart=5.5;
nevol=2000;

%% First phase plane (lifespan=4.5)

subplot(4,3,[7;10])
bJ=1/4.5;
bA=1/4.5;

% Calculate fitness gradients:
[fitgradHval,fitgradPval,tolJvalvec,alphavalvec,~]=JL_fitness_gradients(SSres,tolJmin,alphamin,tolJmax,alphamax,a0,g,q,c1,c2,rJ,rA,tA,f,beta0,bJ,bA,gamma,initvec,t_max);

% Scale fitness gradients according to relative mutation rates:
fitgradHval=fitgradHval*hostmutationprob/(1-hostmutationprob);
tolJvec=tolJvalvec;
alphavec=alphavalvec/alphamax;

% Set up vectors to use later:
xvec=zeros(length(tolJvec)*length(alphavec),1);
yvec=zeros(length(tolJvec)*length(alphavec),1);
plotter1=zeros(length(tolJvec)*length(alphavec),1);
plotter2=zeros(length(tolJvec)*length(alphavec),1);
xcurveJ=zeros(length(tolJvec)*length(alphavec),1);
ycurveJ=zeros(length(alphavec)*length(tolJvec),1);
xcurveA=zeros(length(tolJvec)*length(alphavec),1);
ycurveA=zeros(length(alphavec)*length(tolJvec),1);

% Determine where the fitness gradients change sign:
for i=1:length(tolJvec)
    for j=1:length(alphavec)
        xvec(i+length(tolJvec)*(j-1))=tolJvec(i);
        yvec(i+length(tolJvec)*(j-1))=alphavec(j);
        plotter1(i+length(tolJvec)*(j-1))=fitgradHval(i,j);
        plotter2(i+length(tolJvec)*(j-1))=fitgradPval(i,j);
        
        if j<length(alphavec)
            if fitgradHval(i,j)>0 && fitgradHval(i,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j>1
            if fitgradHval(i,j)>0 && fitgradHval(i,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        
        if j>1 && i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j<length(alphavec) && i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j<length(alphavec) && i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j>1 && i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        
    end
end
xcurveJ=nonzeros(xcurveJ);
ycurveJ=nonzeros(ycurveJ);
xcurveA=nonzeros(xcurveA);
ycurveA=nonzeros(ycurveA);

% Scale the fitness gradient direction arrows so that they are all the same
% length:
for i=1:length(tolJvec)*length(alphavec)
    scaler=sqrt(plotter1(i)^2 + plotter2(i)^2);
    if scaler~=0
        plotter1(i)=plotter1(i)/scaler;
        plotter2(i)=plotter2(i)/scaler;
    end  
end

% Reduce the number of fitness gradient arrows plotted:
deletevec=zeros(length(tolJvec)*length(alphavec),1);
for i=1:length(tolJvec)
    for j=1:length(alphavec)
        if floor((i-(SSres/20))/(SSres/10))~=(i-(SSres/20))/(SSres/10)
            deletevec(i+length(tolJvec)*(j-1))=i+length(tolJvec)*(j-1);
        end
        if floor((j-(SSres/20))/(SSres/10))~=(j-(SSres/20))/(SSres/10)
            deletevec(i+length(tolJvec)*(j-1))=i+length(tolJvec)*(j-1);
        end
    end
end
deletevec=nonzeros(deletevec);
xvec(deletevec)=[];
yvec(deletevec)=[];
plotter1(deletevec)=[];
plotter2(deletevec)=[];

% Create the plot:
qplot=quiver(xvec,yvec,plotter1,plotter2,'color','k');
qplot.AutoScaleFactor=0.2;
xlim([0,1])
ylim([0,1])
set(gca,'xtick',0:0.2:1,'xticklabel',tolJvec(1):(tolJvec(end)-tolJvec(1))/5:tolJvec(end));
set(gca,'ytick',0:0.2:1,'yticklabel',alphamax*(alphavec(1):(alphavec(end)-alphavec(1))/5:alphavec(end)));
hold on
plot(xcurveJ/SSres,ycurveJ/SSres,'o','color',black)
hold on
plot(xcurveA/SSres,ycurveA/SSres,'x','color',grey)
ax=gca;
ax.FontSize=14;
xlabel('Juvenile tolerance, $t_J$','interpreter','latex','fontsize',16)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',16)
axis square
title("lifespan 4.5",'interpreter','latex')
text(0.07,0.93,"C",'fontsize',18)

% Add a trajectory
rng(6)
[trajH,trajP]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,tolJmin,alphamin,tolJmax,alphamax,nevol,hostmutationprob);
trajP=trajP/alphamax;
hold on
% Smooth the data and plot the trajectory:
smooth_trajH=smoothdata(trajH,'movmean',100);
smooth_trajP=smoothdata(trajP,'movmean',100);
plot(smooth_trajH,smooth_trajP,'linewidth',2,'color',red)

%% Second phase plane (lifespan=6.5)

subplot(4,3,[8;11])
bJ=1/6.5;
bA=1/6.5;

% Calculate fitness gradients:
[fitgradHval,fitgradPval,tolJvalvec,alphavalvec,~]=JL_fitness_gradients(SSres,tolJmin,alphamin,tolJmax,alphamax,a0,g,q,c1,c2,rJ,rA,tA,f,beta0,bJ,bA,gamma,initvec,t_max);

% Scale the fitness gradients according to the relative mutation
% probabilities:
fitgradHval=fitgradHval*hostmutationprob/(1-hostmutationprob);
tolJvec=tolJvalvec;
alphavec=alphavalvec/alphamax;

% Set up vectors to use later:
xvec=zeros(length(tolJvec)*length(alphavec),1);
yvec=zeros(length(tolJvec)*length(alphavec),1);
plotter1=zeros(length(tolJvec)*length(alphavec),1);
plotter2=zeros(length(tolJvec)*length(alphavec),1);
xcurveJ=zeros(length(tolJvec)*length(alphavec),1);
ycurveJ=zeros(length(alphavec)*length(tolJvec),1);
xcurveA=zeros(length(tolJvec)*length(alphavec),1);
ycurveA=zeros(length(alphavec)*length(tolJvec),1);

% Determine where the fitness gradients change sign:
for i=1:length(tolJvec)
    for j=1:length(alphavec)
        xvec(i+length(tolJvec)*(j-1))=tolJvec(i);
        yvec(i+length(tolJvec)*(j-1))=alphavec(j);
        plotter1(i+length(tolJvec)*(j-1))=fitgradHval(i,j);
        plotter2(i+length(tolJvec)*(j-1))=fitgradPval(i,j);
        
        if j<length(alphavec)
            if fitgradHval(i,j)>0 && fitgradHval(i,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j>1
            if fitgradHval(i,j)>0 && fitgradHval(i,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        
        if j>1 && i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j<length(alphavec) && i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j<length(alphavec) && i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j>1 && i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        
    end
end
xcurveJ=nonzeros(xcurveJ);
ycurveJ=nonzeros(ycurveJ);
xcurveA=nonzeros(xcurveA);
ycurveA=nonzeros(ycurveA);

% Scale the fitness gradient direction arrows so that they are all the same
% length:
for i=1:length(tolJvec)*length(alphavec)
    scaler=sqrt(plotter1(i)^2 + plotter2(i)^2);
    if scaler~=0
        plotter1(i)=plotter1(i)/scaler;
        plotter2(i)=plotter2(i)/scaler;
    end  
end

% Reduce the number of fitness gradient arrows plotted:
deletevec=zeros(length(tolJvec)*length(alphavec),1);
for i=1:length(tolJvec)
    for j=1:length(alphavec)
        if floor((i-(SSres/20))/(SSres/10))~=(i-(SSres/20))/(SSres/10)
            deletevec(i+length(tolJvec)*(j-1))=i+length(tolJvec)*(j-1);
        end
        if floor((j-(SSres/20))/(SSres/10))~=(j-(SSres/20))/(SSres/10)
            deletevec(i+length(tolJvec)*(j-1))=i+length(tolJvec)*(j-1);
        end
    end
end
deletevec=nonzeros(deletevec);
xvec(deletevec)=[];
yvec(deletevec)=[];
plotter1(deletevec)=[];
plotter2(deletevec)=[];
        
% Create the plot:
qplot=quiver(xvec,yvec,plotter1,plotter2,'color','k');
qplot.AutoScaleFactor=0.2;
xlim([0,1])
ylim([0,1])
set(gca,'xtick',0:0.2:1,'xticklabel',tolJvec(1):(tolJvec(end)-tolJvec(1))/5:tolJvec(end));
set(gca,'ytick',0:0.2:1,'yticklabel',alphamax*(alphavec(1):(alphavec(end)-alphavec(1))/5:alphavec(end)));
hold on
plot(xcurveJ/SSres,ycurveJ/SSres,'o','color',black)
hold on
plot(xcurveA/SSres,ycurveA/SSres,'x','color',grey)
ax=gca;
ax.FontSize=14;
xlabel('Juvenile tolerance, $t_J$','interpreter','latex','fontsize',16)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',16)
axis square
title("lifespan 6.5",'interpreter','latex')
text(0.07,0.93,"D",'fontsize',18)

% Add a trajectory
rng(4)
[trajH,trajP]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,tolJmin,alphamin,tolJmax,alphamax,nevol,hostmutationprob);
trajP=trajP/alphamax;
hold on
% Smooth the data and plot the trajectory:
smooth_trajH=smoothdata(trajH,'movmean',100);
smooth_trajP=smoothdata(trajP,'movmean',100);
plot(smooth_trajH,smooth_trajP,'linewidth',2,'color',red)

%% Third phase plane (lifespan=8.5)

subplot(4,3,[9;12])
bJ=1/8.5;
bA=1/8.5;

% Calculate the fitness gradients:
[fitgradHval,fitgradPval,tolJvalvec,alphavalvec,~]=JL_fitness_gradients(SSres,tolJmin,alphamin,tolJmax,alphamax,a0,g,q,c1,c2,rJ,rA,tA,f,beta0,bJ,bA,gamma,initvec,t_max);

% Scale the fitness gradients according to the relative mutation
% probabilities:
fitgradHval=fitgradHval*hostmutationprob/(1-hostmutationprob);
tolJvec=tolJvalvec;
alphavec=alphavalvec/alphamax;

% Set up vectors to use later:
xvec=zeros(length(tolJvec)*length(alphavec),1);
yvec=zeros(length(tolJvec)*length(alphavec),1);
plotter1=zeros(length(tolJvec)*length(alphavec),1);
plotter2=zeros(length(tolJvec)*length(alphavec),1);
xcurveJ=zeros(length(tolJvec)*length(alphavec),1);
ycurveJ=zeros(length(alphavec)*length(tolJvec),1);
xcurveA=zeros(length(tolJvec)*length(alphavec),1);
ycurveA=zeros(length(alphavec)*length(tolJvec),1);

% Determine where the fitness gradients change sign:
for i=1:length(tolJvec)
    for j=1:length(alphavec)
        xvec(i+length(tolJvec)*(j-1))=tolJvec(i);
        yvec(i+length(tolJvec)*(j-1))=alphavec(j);
        plotter1(i+length(tolJvec)*(j-1))=fitgradHval(i,j);
        plotter2(i+length(tolJvec)*(j-1))=fitgradPval(i,j);
        
        if j<length(alphavec)
            if fitgradHval(i,j)>0 && fitgradHval(i,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j>1
            if fitgradHval(i,j)>0 && fitgradHval(i,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        
        if j>1 && i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j<length(alphavec) && i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j<length(alphavec) && i>1
            if fitgradHval(i,j)>0 && fitgradHval(i-1,j+1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i-1,j+1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        if j>1 && i<length(tolJvec)
            if fitgradHval(i,j)>0 && fitgradHval(i+1,j-1)<=0
                xcurveJ(i+length(tolJvec)*(j-1))=i;
                ycurveJ(i+length(tolJvec)*(j-1))=j;
            end
            if fitgradPval(i,j)>0 && fitgradPval(i+1,j-1)<=0
                xcurveA(i+length(tolJvec)*(j-1))=i;
                ycurveA(i+length(tolJvec)*(j-1))=j;
            end
        end
        
    end
end
xcurveJ=nonzeros(xcurveJ);
ycurveJ=nonzeros(ycurveJ);
xcurveA=nonzeros(xcurveA);
ycurveA=nonzeros(ycurveA);

% Scale the fitness gradient arrows so that they are all the same length:
for i=1:length(tolJvec)*length(alphavec)
    scaler=sqrt(plotter1(i)^2 + plotter2(i)^2);
    if scaler~=0
        plotter1(i)=plotter1(i)/scaler;
        plotter2(i)=plotter2(i)/scaler;
    end  
end

% Reduce the number of fitness gradient arrows plotted:
deletevec=zeros(length(tolJvec)*length(alphavec),1);
for i=1:length(tolJvec)
    for j=1:length(alphavec)
        if floor((i-(SSres/20))/(SSres/10))~=(i-(SSres/20))/(SSres/10)
            deletevec(i+length(tolJvec)*(j-1))=i+length(tolJvec)*(j-1);
        end
        if floor((j-(SSres/20))/(SSres/10))~=(j-(SSres/20))/(SSres/10)
            deletevec(i+length(tolJvec)*(j-1))=i+length(tolJvec)*(j-1);
        end
    end
end
deletevec=nonzeros(deletevec);
xvec(deletevec)=[];
yvec(deletevec)=[];
plotter1(deletevec)=[];
plotter2(deletevec)=[];

% Create the plot:
qplot=quiver(xvec,yvec,plotter1,plotter2,'color','k');
qplot.AutoScaleFactor=0.2;
xlim([0,1])
ylim([0,1])
set(gca,'xtick',0:0.2:1,'xticklabel',tolJvec(1):(tolJvec(end)-tolJvec(1))/5:tolJvec(end));
set(gca,'ytick',0:0.2:1,'yticklabel',alphamax*(alphavec(1):(alphavec(end)-alphavec(1))/5:alphavec(end)));
hold on
plot(xcurveJ/SSres,ycurveJ/SSres,'o','color',black)
hold on
plot(xcurveA/SSres,ycurveA/SSres,'x','color',grey)
ax=gca;
ax.FontSize=14;
xlabel('Juvenile tolerance, $t_J$','interpreter','latex','fontsize',16)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',16)
axis square
title("lifespan 8.5",'interpreter','latex')
text(0.07,0.93,"E",'fontsize',18)

% Add a trajectory
rng(9)
nevol=20000;
alpha_max=25;
[trajH,trajP]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,tolJmin,alphamin,tolJmax,alpha_max,nevol,hostmutationprob);
trajP=trajP/alphamax;
hold on
% Smooth the data and plot the trajectory:
smooth_trajH=smoothdata(trajH,'movmean',100);
smooth_trajP=smoothdata(trajP,'movmean',100);
smooth_trajP(isnan(smooth_trajH))=[];
smooth_trajH(isnan(smooth_trajH))=[];
smooth_trajP(642:1162)=[];
smooth_trajH(642:1162)=[];
plot(smooth_trajH,smooth_trajP,'linewidth',2,'color',red)
