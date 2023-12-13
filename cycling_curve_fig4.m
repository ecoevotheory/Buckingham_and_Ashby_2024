% This code uses simulations to determine whether cycling occurs for
% different values of the host lifespan and then creates a plot showing the
% maximum and minimum extent of the cycling trait for each value of the
% host lifespan. 

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
    init_pop = [0.1,0.1,0.1,0.1];
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
    [tolJ_start,alpha_start,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,TOLJ,ALPHA,DISPREV,N] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tolJmin,tolJmax,tolJ_start,alphamin,alphamax,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
    
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

% Colours for plotting:
red=1/255*[215,48,39];
orange=1/255*[253,174,97];
blue=1/255*[69,117,180];

% Smooth parts of the curves:
smooth_traj1=smoothdata(cycle_max_vec,'movmean',5);
smooth_traj2=smoothdata(cycle_min_vec,'movmean',5);
smooth_traj1=[cycle_max_vec(1:4);smooth_traj1(5:16);cycle_max_vec(17:end)];
smooth_traj2=[cycle_min_vec(1:4);smooth_traj2(5:16);cycle_min_vec(17:end)];
cycle_both_vec(5)=NaN;

% Make the plot for virulence:
subplot(2,1,2)
plot(lifespanvec,smooth_traj1,'color',red,'linewidth',3)
hold on
plot(lifespanvec,smooth_traj2,'color',orange,'linewidth',3)
hold on
plot(lifespanvec,cycle_both_vec,'color',blue,'linewidth',3)
xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Virulence, $\alpha$','interpreter','latex')
ylim([0,10])
xlim([4,9])
text(8.85,8.6,"B",'fontsize',18)
set(gca,'xtick',[4,5,6,7,8,9],'fontsize',14)
set(gca,'ytick',[0,5,10])
pbaspect([2,1,1])

% Smooth parts of the curves:
smooth_traj3=smoothdata(tol_max_vec,'movmean',5);
smooth_traj4=smoothdata(tol_min_vec,'movmean',5);
smooth_traj3=[tol_max_vec(1:4);smooth_traj3(5:16);tol_max_vec(17:end)];
smooth_traj4=[tol_min_vec(1:4);smooth_traj4(5:16);tol_min_vec(17:end)];
tol_both_vec(5)=NaN;

% Make the plot for tolerance:
subplot(2,1,1)
plot(lifespanvec,smooth_traj3,'color',red,'linewidth',3)
hold on
plot(lifespanvec,smooth_traj4,'color',orange,'linewidth',3)
hold on
plot(lifespanvec,tol_both_vec,'color',blue,'linewidth',3)
xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Juvenile tolerance, $\tau_J$','interpreter','latex')
ylim([0.5,1])
xlim([4,9])
text(8.85,0.92,"A",'fontsize',18)
set(gca,'ytick',[0.5,0.75,1])
set(gca,'xtick',[4,5,6,7,8,9],'fontsize',14)
pbaspect([2,1,1])
