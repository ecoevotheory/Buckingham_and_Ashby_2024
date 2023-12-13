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
cycle_both_vec(5)=NaN;

% Plot the virulence:
subplot(4,3,[4;6])
plot(lifespanvec,smooth_traj1,'color',red,'linewidth',2)
hold on
plot(lifespanvec,smooth_traj2,'color',orange,'linewidth',2)
hold on
plot(lifespanvec,cycle_both_vec,'color',blue,'linewidth',2)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',14)
ylim([0,10])
xlim([4,9])
text(8.85,8.6,"B",'fontsize',18)
set(gca,'xtick',[4,5,6,7,8,9],'fontsize',14)
xlabel('Lifespan, $1/b$','interpreter','latex','fontsize',14)

% Smooth parts of the curves:
smooth_traj3=smoothdata(tol_max_vec,'movmean',5);
smooth_traj4=smoothdata(tol_min_vec,'movmean',5);
smooth_traj3=[tol_max_vec(1:4);smooth_traj3(5:16);tol_max_vec(17:end)];
smooth_traj4=[tol_min_vec(1:4);smooth_traj4(5:16);tol_min_vec(17:end)];
tol_both_vec(5)=NaN;

% Plot the tolerance:
subplot(4,3,[1;3])
plot(lifespanvec,smooth_traj3,'color',red,'linewidth',2)
hold on
plot(lifespanvec,smooth_traj4,'color',orange,'linewidth',2)
hold on
plot(lifespanvec,tol_both_vec,'color',blue,'linewidth',2)
xlabel('Lifespan, $1/b$','interpreter','latex','fontsize',14)
ylabel('Tolerance, $t_J$','interpreter','latex','fontsize',14)
ylim([0.5,1])
xlim([4,9])
text(8.85,0.92,"A",'fontsize',18)
set(gca,'ytick',[0.5,0.75,1])
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
rng(2)
[trajH,trajP]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,tolJmin,alphamin,tolJmax,alphamax,nevol,hostmutationprob);
trajP=trajP/alphamax;
hold on
% Smooth the data and plot the trajectory:
smooth_trajH=smoothdata(trajH,'movmean',100);
smooth_trajP=smoothdata(trajP,'movmean',100);
plot(smooth_trajH,smooth_trajP,'linewidth',2,'color',red)
