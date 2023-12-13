% This code draws two phase planes with example trajectories.

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

initvec=[0.1,0.1,0.01,0.01];
orig_tmax=100;
maxsingstrats=1000;
SSres=500;
startH=0;
startP=0;
finH=1;
finP=20;

% Find fitness gradients for different values of host and parasite traits:
[fitgradHval,fitgradPval,tolJvalvec,alphavalvec,R0counter]=JL_fitness_gradients(SSres,startH,startP,finH,finP,a0,g,q,c1,c2,rJ,rA,tA,f,beta0,bJ,bA,gamma,initvec,orig_tmax);

%% Find and plot nullclines and fitness gradient directions

% Set up vectors:
tolJvec=tolJvalvec;
alphavec=alphavalvec/finP;

xvec=zeros(length(tolJvec)*length(alphavec),1);
yvec=zeros(length(tolJvec)*length(alphavec),1);
plotter1=zeros(length(tolJvec)*length(alphavec),1);
plotter2=zeros(length(tolJvec)*length(alphavec),1);
xcurveJ=zeros(length(tolJvec)*length(alphavec),1);
ycurveJ=zeros(length(alphavec)*length(tolJvec),1);
xcurveA=zeros(length(tolJvec)*length(alphavec),1);
ycurveA=zeros(length(alphavec)*length(tolJvec),1);

% Find where fitness gradients change sign:
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

% Scale fitness gradient arrows so that they are all the same length:
for i=1:length(tolJvec)*length(alphavec)
    scaler=sqrt(plotter1(i)^2 + plotter2(i)^2);
    if scaler~=0
        plotter1(i)=plotter1(i)/scaler;
        plotter2(i)=plotter2(i)/scaler;
    end  
end

% Reduce the number of fitness gradient arrows on the plot:
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
        
% Colours for plotting:
red=1/255*[215,48,39];
blue=1/255*[69,117,180];
orange=1/255*[253,174,97];
black=1/255*[0,0,0];
grey=1/255*[166,166,166];

% Create the plot:
subplot(1,2,1)
qplot=quiver(xvec,yvec,plotter1,plotter2,'color','k');
qplot.AutoScaleFactor=0.2;
xlim([0,1])
ylim([0,1])
set(gca,'xtick',0:0.2:1,'xticklabel',tolJvec(1):(tolJvec(end)-tolJvec(1))/5:tolJvec(end));
set(gca,'ytick',0:0.2:1,'yticklabel',finP*(alphavec(1):(alphavec(end)-alphavec(1))/5:alphavec(end)));
hold on
plot(xcurveJ/SSres,ycurveJ/SSres,'o','color',black)
hold on
plot(xcurveA/SSres,ycurveA/SSres,'x','color',grey)
ax=gca;
ax.FontSize=14;
xlabel('Juvenile tolerance, $t_J$','interpreter','latex','fontsize',16)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',16)
title('A','fontsize',18)
axis square

subplot(1,2,2)
qplot=quiver(xvec,yvec,plotter1,plotter2,'color','k');
qplot.AutoScaleFactor=0.2;
xlim([0,1])
ylim([0,1])
set(gca,'xtick',0:0.2:1,'xticklabel',tolJvec(1):(tolJvec(end)-tolJvec(1))/5:tolJvec(end));
set(gca,'ytick',0:0.2:1,'yticklabel',finP*(alphavec(1):(alphavec(end)-alphavec(1))/5:alphavec(end)));
hold on
h1=plot(xcurveJ/SSres,ycurveJ/SSres,'o','color',black);
hold on
h2=plot(xcurveA/SSres,ycurveA/SSres,'x','color',grey);
ax=gca;
ax.FontSize=14;
xlabel('Juvenile tolerance, $t_J$','interpreter','latex','fontsize',16)
ylabel('Virulence, $\alpha$','interpreter','latex','fontsize',16)
title('B','fontsize',18)
axis square

%% Add trajectories

% Pick a starting point and run a simulation:
subplot(1,2,1)
rng(3)
tolstart=0.9;
alphastart=10;
nevol=5000;
hostmutationprob=0.5;
[trajH,trajP]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,startH,startP,finH,finP,nevol,hostmutationprob);
trajP=trajP/finP;
hold on
% Smooth the data and plot the trajectory:
smooth_trajH=smoothdata(trajH,'movmean',100);
smooth_trajP=smoothdata(trajP,'movmean',100);
h3=plot(smooth_trajH,smooth_trajP,'linewidth',2,'color',blue);

%% Pick another starting point and run a simulation:
subplot(1,2,1)
rng(1)
tolstart=0.9;
alphastart=10;
nevol=2000;
hostmutationprob=1/11;
[trajH2,trajP2]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,startH,startP,finH,finP,nevol,hostmutationprob);
trajP2=trajP2/finP;
deletevec=zeros(nevol,1);
for i=1:nevol
    if isnan(trajH2(i)) || isnan(trajP2(i))
        deletevec(i)=i;
    end
end
deletevec=nonzeros(deletevec);
trajH2(deletevec)=[];
trajP2(deletevec)=[];
hold on
% Smooth the data and plot the trajectory:
smooth_trajH2=smoothdata(trajH2,'movmean',100);
smooth_trajP2=smoothdata(trajP2,'movmean',100);
h4=plot(smooth_trajH2,smooth_trajP2,'linewidth',2,'color',red);

%% Pick another starting point and run a simulation:
subplot(1,2,2)
rng(3)
tolstart=0.9;
alphastart=3.7;
nevol=1200;
hostmutationprob=0.5;
[trajH3,trajP3]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,startH,startP,finH,finP,nevol,hostmutationprob);
trajP3=trajP3/finP;
deletevec=zeros(nevol,1);
for i=1:nevol
    if isnan(trajH3(i)) || isnan(trajP3(i))
        deletevec(i)=i;
    end
end
deletevec=nonzeros(deletevec);
trajH3(deletevec)=[];
trajP3(deletevec)=[];
hold on
% Smooth the data and plot the trajectory:
smooth_trajH3=smoothdata(trajH3,'movmean',100);
smooth_trajP3=smoothdata(trajP3,'movmean',100);
plot(smooth_trajH3,smooth_trajP3,'linewidth',2,'color',blue)

%% Pick another starting point and run a simulation:
subplot(1,2,2)
rng(2)
tolstart=0.9;
alphastart=3.7;
nevol=1800;
hostmutationprob=1/11;
[trajH4,trajP4]=JL_sim_traj_function(tolstart,alphastart,q,f,bJ,g,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,startH,startP,finH,finP,nevol,hostmutationprob);
trajP4=trajP4/finP;
deletevec=zeros(nevol,1);
for i=1:nevol
    if isnan(trajH4(i)) || isnan(trajP4(i))
        deletevec(i)=i;
    end
end
deletevec=nonzeros(deletevec);
trajH4(deletevec)=[];
trajP4(deletevec)=[];
hold on
% Smooth the data and plot the trajectory:
smooth_trajH4=smoothdata(trajH4,'movmean',100);
smooth_trajP4=smoothdata(trajP4,'movmean',100);
plot(smooth_trajH4,smooth_trajP4,'linewidth',2,'color',red)

%% Indicate the starting-points of the trajectories:
subplot(1,2,1)
hold on
plot(smooth_trajH(1)*0.5+smooth_trajH2(1)*0.5,smooth_trajP(1)*0.5+smooth_trajP2(1)*0.5,'o','markersize',12,'color',orange,'linewidth',4)
subplot(1,2,2)
hold on
plot(smooth_trajH3(1),smooth_trajP3(1),'o','markersize',12,'color',orange,'linewidth',4)
