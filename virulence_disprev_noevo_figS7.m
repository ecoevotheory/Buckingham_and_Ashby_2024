% This code shows the effect of increasing pathogen virulence on the level
% of disease prevalence (when tolerance is held evolutionarily static). 

% Define parameter values:
finP=10;
numberofpoints=50;
vir_vec=linspace(finP/numberofpoints,finP,numberofpoints);
tol_vec=[0.2,0.5,0.7];

rJ=0;
rA=0;
g=1;
a0=5;
c1=0.25;
c2=3;
q=1;
beta0=10;
bJ=1/5;
bA=1/5;
f=1;
gamma=0;
initvec=[0.1,0.1,0.1,0.1];
orig_tmax=100;
eps=0.001;
index=1;

% Determine the disease prevalence for each combination of tolerance and
% virulence:
LL_disprev=NaN(length(vir_vec),length(tol_vec));
JL_disprev=NaN(length(vir_vec),length(tol_vec));
for j=1:length(tol_vec)
    tol=tol_vec(j);
    a=a0*(1-c1*(1-exp(c2*tol))/(1-exp(c2)));
    
    for i=1:length(vir_vec)
        alpha=vir_vec(i);
        
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tol,tol,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
        
        popdens=SJval+SAval+IAval+IJval;
        
        if popdens>eps
            LL_disprev(i,j)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        end
        
        [SJval,SAval,IJval,IAval,~]=endemic_equilibrium_function(tol,0,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
        
        popdens=SJval+SAval+IAval+IJval;
        
        if popdens>eps
            JL_disprev(i,j)=(IJval+IAval)/(SJval+SAval+IJval+IAval);
        end
        
    end
end

% Find the proportional reduction in disease prevalence:
LL_disprev_prop_reduc=NaN(length(vir_vec),length(tol_vec));
JL_disprev_prop_reduc=NaN(length(vir_vec),length(tol_vec));
for j=1:length(tol_vec)
    for i=1:length(vir_vec)-1
        LL_disprev_prop_reduc(i+1,j)=(LL_disprev(index,j)-LL_disprev(i+1,j))/LL_disprev(index,j);
        JL_disprev_prop_reduc(i+1,j)=(JL_disprev(index,j)-JL_disprev(i+1,j))/JL_disprev(index,j);
    end
end
LL_disprev_prop_reduc(1,:)=zeros(1,length(tol_vec));
JL_disprev_prop_reduc(1,:)=zeros(1,length(tol_vec));

% Create the plot

blue=1/255*[69,117,180];
red=1/255*[215,48,39];

hold on
plot(vir_vec,LL_disprev_prop_reduc(:,1),'color',blue,'linestyle','-','linewidth',2)
plot(vir_vec,LL_disprev_prop_reduc(:,2),'color',blue,'linestyle','--','linewidth',2)
plot(vir_vec,LL_disprev_prop_reduc(:,3),'color',blue,'linestyle',':','linewidth',2)
plot(vir_vec,JL_disprev_prop_reduc(:,1),'color',red,'linestyle','-','linewidth',2)
plot(vir_vec,JL_disprev_prop_reduc(:,2),'color',red,'linestyle','--','linewidth',2)
plot(vir_vec,JL_disprev_prop_reduc(:,3),'color',red,'linestyle',':','linewidth',2)

xlabel('Pathogen mortality virulence, $\alpha$','interpreter','latex','fontsize',18)
ylabel('Proportional reduction in disease prevalence','interpreter','latex','fontsize',18)
xlim([0,finP])
ylim([min(LL_disprev_prop_reduc(:,3)),1])
set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1])
