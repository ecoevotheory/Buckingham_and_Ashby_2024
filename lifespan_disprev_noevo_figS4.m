% This codes plots the effect of lifespan on disease prevalence when both
% tolerance and virulence are held evolutionary static. 

% Use these values for lifespan, tolerance and virulence:
lifespanvec=[1,1.3,1.325,1.33,1.34,1.35,1.4,1.45,1.5,1.75,2,2.25,2.5,3,4,5,6,7,8,9,10];
LL_tol_vec=[0.2,0.2,0.2,0.7,0.7,0.7];
LL_vir_vec=[1,3,5,1,3,5];
JL_tol_vec=[0,0,0,0.2,0.2,0.2];
JL_vir_vec=[1,3,5,1,3,5];

% Define parameter values:
rJ=0;
rA=0;
g=1;
a0=1;
c1=1;
c2=4;
q=1;
beta0=10;
f=1;
gamma=1;
initvec=[0.1,0.1,0.1,0.1];
orig_tmax=100;
eps=0.001;

% Calculate the disease prevalence for each combination of tolerance and
% virulence in the lifelong tolerance scenario:
disprev_LL=NaN(length(lifespanvec),6);
for j=1:6
    tol=LL_tol_vec(j);
    alpha=LL_vir_vec(j);
    a=a0*(1-c1*(1-exp(c2*tol))/(1-exp(c2)));
    for i=1:length(lifespanvec)
        
        lifespan=lifespanvec(i);
        bJ=1/lifespan;
        bA=1/lifespan;
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tol,tol,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
        
        pop_dens=SJval+SAval+IJval+IAval;
        
        if pop_dens>eps
            disprev_LL(i,j)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        end
        
    end
end

% Calculate the disease prevalence for each combination of tolerance and
% virulence in the juvenile tolerance scenario:
disprev_JL=NaN(length(lifespanvec),6);
for j=1:6
    tol=JL_tol_vec(j);
    alpha=JL_vir_vec(j);
    a=a0*(1-c1*(1-exp(c2*tol))/(1-exp(c2)));
    for i=1:length(lifespanvec)
        
        lifespan=lifespanvec(i);
        bJ=1/lifespan;
        bA=1/lifespan;
        
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tol,0,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
        
        pop_dens=SJval+SAval+IJval+IAval;
        
        if pop_dens>eps
            disprev_JL(i,j)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        end
        
        
    end
end

% Create the plot

blue=1/255*[69,117,180];
red=1/255*[215,48,39];
orange=1/255*[253,174,97];

subplot(1,2,1)
hold on
plot(lifespanvec,disprev_LL(:,1),'color',blue,'linestyle','-','linewidth',2)
plot(lifespanvec,disprev_LL(:,2),'color',orange,'linestyle','-','linewidth',2)
plot(lifespanvec,disprev_LL(:,3),'color',red,'linestyle','-','linewidth',2)
plot(lifespanvec,disprev_LL(:,4),'color',blue,'linestyle','--','linewidth',2)
plot(lifespanvec,disprev_LL(:,5),'color',orange,'linestyle','--','linewidth',2)
plot(lifespanvec,disprev_LL(:,6),'color',red,'linestyle','--','linewidth',2)
xlim([0,10])
ylim([0,0.5])
xlabel('Lifespan, $1/b$','interpreter','latex','fontsize',18)
ylabel('Disease prevalence','interpreter','latex','fontsize',18)
box on
set(gca,'ytick',[0,0.1,0.2,0.3,0.4,0.5],'fontsize',16)
text(0.6,0.47,'A','fontsize',30)

subplot(1,2,2)
hold on
plot(lifespanvec,disprev_JL(:,1),'color',blue,'linestyle','-','linewidth',2)
plot(lifespanvec,disprev_JL(:,2),'color',orange,'linestyle','-','linewidth',2)
plot(lifespanvec,disprev_JL(:,3),'color',red,'linestyle','-','linewidth',2)
plot(lifespanvec,disprev_JL(:,4),'color',blue,'linestyle','--','linewidth',2)
plot(lifespanvec,disprev_JL(:,5),'color',orange,'linestyle','--','linewidth',2)
plot(lifespanvec,disprev_JL(:,6),'color',red,'linestyle','--','linewidth',2)
xlim([0,10])
ylim([0,0.5])
xlabel('Lifespan, $1/b$','interpreter','latex','fontsize',18)
ylabel('Disease prevalence','interpreter','latex','fontsize',18)
box on
set(gca,'ytick',[0,0.1,0.2,0.3,0.4,0.5],'fontsize',16)
text(0.6,0.47,'B','fontsize',30)

% colour indicates level of virulence (blue is low, orange is medium and
% red is high)
% line style indicates level of tolerance (solid is low and dashed is high)
