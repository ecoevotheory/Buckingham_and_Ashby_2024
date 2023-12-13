function deriv=JL_find_parasite_singstrat_derivative(q,f,g,bJ,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,tJval,initvec,orig_tmax)

% This function finds the derivative of the parasite singular strategy
% (written as a function of the host trait in the case where the host
% evolves much more slowly than the pathogen). 

eps=0.01;

% Find the value of singular value of the pathogen trait when the host 
% trait is just below the value at which we are evaluating the derivative:
tJ=tJval-eps;
a=a0*(1-(c1*(1-exp(c2*tJ)))/(1-exp(c2)));
parasite_fitgrad_sign=NaN(10000,1);
findervec=NaN(9999,1);
for i=1:10000
    alpha=100*i/10000;
    [SJ,SA,~,~,~]=endemic_equilibrium_function(tJ,tA,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);           
    parasite_fitgrad_sign(i)=0.5*(1-rJ)*SJ*(bA+alpha*(1-tA)+g+gamma)*(bA+alpha*(1-tA)+gamma)*(bJ+g+gamma-alpha*(1-tJ))+0.5*(1-rA)*SA*(bJ+g+gamma+alpha*(1-tJ))^2*(bA-alpha*(1-tA)+gamma)-(1-tA)*(1-rJ)*SJ*alpha*g*(bJ+g+alpha*(1-tJ)+gamma);
    if i>1
        findervec(i)=parasite_fitgrad_sign(i)*parasite_fitgrad_sign(i-1);
    end
end         
alpha_below=100*(find(findervec<0))/10000;

% Find the value of singular value of the pathogen trait when the host 
% trait is just above the value at which we are evaluating the derivative:
tJ=tJval+eps;
a=a0*(1-(c1*(1-exp(c2*tJ)))/(1-exp(c2)));
parasite_fitgrad_sign=NaN(10000,1);
findervec=NaN(9999,1);
for i=1:10000
    alpha=100*i/10000;
    [SJ,SA,~,~,~]=endemic_equilibrium_function(tJ,tA,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax);
    parasite_fitgrad_sign(i)=0.5*(1-rJ)*SJ*(bA+alpha*(1-tA)+g+gamma)*(bA+alpha*(1-tA)+gamma)*(bJ+g+gamma-alpha*(1-tJ))+0.5*(1-rA)*SA*(bJ+g+gamma+alpha*(1-tJ))^2*(bA-alpha*(1-tA)+gamma)-(1-tA)*(1-rJ)*SJ*alpha*g*(bJ+g+alpha*(1-tJ)+gamma);
    if i>1
        findervec(i)=parasite_fitgrad_sign(i)*parasite_fitgrad_sign(i-1);
    end
end         
alpha_above=100*(find(findervec<0))/10000;

% The derivative is given by:
if length(alpha_above)~=length(alpha_below)
    deriv=[];
else
    deriv=(alpha_above-alpha_below)/(2*eps);
end
if length(deriv)>1
    deriv=[];
end

end