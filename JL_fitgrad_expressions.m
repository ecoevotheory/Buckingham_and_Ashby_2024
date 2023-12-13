% This code determines analytical expressions for the fitness gradient and
% the second derivatives of the invasion fitness. These expressions are 
% used directly throughout the code for this model.

% Define variables and trade-offs:
syms SJ SA IJ IA tJ resA alpha tJm alpham tA rJ rA g gamma a0 f c1 c2 beta0 bJ bA tJ alpha q
a=a0*(1-(c1*(1-exp(c2*tJ))/(1-exp(c2))));
am=subs(a,tJ,tJm);

% Set up parameters and functions:
lambdaJstar=beta0*(1-rJ)*(sqrt(alpha)*IJ+sqrt(alpha)*IA);
lambdaAstar=beta0*(1-rA)*(sqrt(alpha)*IJ+sqrt(alpha)*IA);
Am=g*((bJ+g+alpha*(1-tJm)+gamma)*(bA+alpha*(1-tA)+gamma)+f*(bJ+g+alpha*(1-tJm)+gamma)*lambdaAstar+f*lambdaJstar*(bA+lambdaAstar))+gamma*lambdaJstar;
Bm=(bJ+g+lambdaJstar)*((bA+alpha*(1-tA)+gamma)*(bJ+g+alpha*(1-tJm)+gamma)*(bA+lambdaAstar)-gamma*lambdaAstar*(bJ+g+alpha*(1-tJm)+gamma))-gamma*((bA+lambdaAstar)*lambdaJstar*(bA+alpha*(1-tA)+gamma)-gamma*lambdaJstar*lambdaAstar);
% Define the invasion fitness:
wH=((am*(1-q*(SJ+SA+IJ+IA))*Am)/Bm) -1;
wP=((beta0*(1-rJ)*SJ*(sqrt(alpham)*(bA+gamma+alpham*(1-tA))+sqrt(alpham)*g)+beta0*(1-rA)*SA*sqrt(alpham)*(bJ+gamma+g+alpham*(1-tJ)))/((bJ+gamma+alpham*(1-tJ)+g)*(bA+gamma+alpham*(1-tA))))-1;
% Determine the fitness gradient:
fitgradHworking=diff(wH,tJm);
fitgradPworking=diff(wP,alpham);
fitgradHfinal=subs(fitgradHworking,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
fitgradPfinal=subs(fitgradPworking,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
% Determine the second derivatives with respect to the mutant traits only:
deriv2Hworking=diff(fitgradHworking,tJm);
deriv2H=subs(deriv2Hworking,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
deriv2Pworking=diff(fitgradPworking,alpham);
deriv2P=subs(deriv2Pworking,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
% Determine the cross-derivatives at (tJ, alpha) and at points nearby:
syms SJupH SAupH IJupH IAupH SJdownH SAdownH IJdownH IAdownH eps SJupP SAupP IJupP IAupP SJdownP SAdownP IJdownP IAdownP
cderivH1working=(subs(fitgradHworking,[tJ,SJ,SA,IJ,IA],[tJ+eps,SJupH,SAupH,IJupH,IAupH])-subs(fitgradHworking,[tJ,SJ,SA,IJ,IA],[tJ-eps,SJdownH,SAdownH,IJdownH,IAdownH]))/(2*eps);
cderivH2working=(subs(fitgradHworking,[alpha,SJ,SA,IJ,IA],[alpha+eps,SJupP,SAupP,IJupP,IAupP])-subs(fitgradHworking,[alpha,SJ,SA,IJ,IA],[alpha-eps,SJdownP,SAdownP,IJdownP,IAdownP]))/(2*eps);
cderivH1=subs(cderivH1working,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
cderivH2=subs(cderivH2working,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
cderivP1working=(subs(fitgradPworking,[alpha,SJ,SA,IJ,IA],[alpha+eps,SJupP,SAupP,IJupP,IAupP])-subs(fitgradPworking,[alpha,SJ,SA,IJ,IA],[alpha-eps,SJdownP,SAdownP,IJdownP,IAdownP]))/(2*eps);
cderivP2working=(subs(fitgradPworking,[tJ,SJ,SA,IJ,IA],[tJ+eps,SJupH,SAupH,IJupH,IAupH])-subs(fitgradPworking,[tJ,SJ,SA,IJ,IA],[tJ-eps,SJdownH,SAdownH,IJdownH,IAdownH]))/(2*eps);
cderivP1=subs(cderivP1working,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);
cderivP2=subs(cderivP2working,[tJm,tJ,alpham,alpha],[tJ,tJ,alpha,alpha]);

% Display the output expressions:
disp("fitgradH=")
disp(fitgradHfinal)
disp("fitgradP=")
disp(fitgradPfinal)
disp("deriv2H=")
disp(deriv2H)
disp("deriv2P=")
disp(deriv2P)
disp("cderivH1=")
disp(cderivH1)
disp("cderivH2=")
disp(cderivH2)
disp("cderivP1=")
disp(cderivP1)
disp("cderivP2=")
disp(cderivP2)
