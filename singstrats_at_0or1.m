function [tJss,alphass]=singstrats_at_0or1(fitgradHval,fitgradPval,tJvalvec,alphavalvec,SSres,maxsingstrats,R0counter,finP)

% This function determines whether there are singular strategies when the
% evolving traits are at their maximum or minimum values (and hence where 
% there may not be sign changes in both fitness gradients).

% Set up parameters and vectors to use later:
finished=0;
tJss=-10+zeros(1,maxsingstrats);
alphass=-10+zeros(1,maxsingstrats);
foundsofar=1;

for j=1:SSres-1
    % Determine fitness gradients for each value of virulence when tolerance is zero:
    if fitgradHval(1,j)<0 && fitgradPval(1,j)>0 && fitgradHval(1,j+1)<0 && fitgradPval(1,j+1)<0 && tJvalvec(1)==0
        newSS=(alphavalvec(j)+alphavalvec(j+1))/2;
        tJss(foundsofar)=0;
        alphass(foundsofar)=newSS;
        foundsofar=foundsofar+1;
        finished=1;
    % Determine fitness gradients for each value of virulence when tolerance is one:
    elseif fitgradHval(end,j)>0 && fitgradPval(end,j)>0 && fitgradHval(end,j+1)>0 && fitgradPval(end,j+1)<0 && tJvalvec(end)==1
        newSS=(alphavalvec(j)+alphavalvec(j+1))/2;
        tJss(foundsofar)=1;
        alphass(foundsofar)=newSS;
        foundsofar=foundsofar+1;
        finished=1;
    % Determine fitness gradients for each value of tolerance when virulence is zero:
    elseif fitgradHval(j,1)>0 && fitgradHval(j+1,1)<0 && fitgradPval(j,1)<0 && fitgradPval(j+1,1)<0 && alphavalvec(1)==0
        newSS=(tJvalvec(j)+tJvalvec(j+1))/2;
        tJss(foundsofar)=newSS;
        alphass(foundsofar)=0;
        foundsofar=foundsofar+1;
        finished=1;
    % Determine fitness gradients for each value of tolerance when virulence is arbitrarily large:
    elseif fitgradHval(j,end)>0 && fitgradHval(j+1,end)<0 && fitgradPval(j,end)>0 && fitgradPval(j+1,end)>0 && alphavalvec(end)==finP
        newSS=(tJvalvec(j)+tJvalvec(j+1))/2;
        tJss(foundsofar)=newSS;
        alphass(foundsofar)=finP;
        foundsofar=foundsofar+1;
        finished=1;
    end
end

% Other cases:
if finished==0
    if sum(isnan(fitgradHval),'all')==SSres^2 % hosts are never viable
        tJss=-10+zeros(1,maxsingstrats);
        alphass=-10+zeros(1,maxsingstrats);
    elseif R0counter==SSres^2 % The disease is never viable
        tJss=[0,-10+zeros(1,maxsingstrats-1)];
        alphass=[-10,-10+zeros(1,maxsingstrats-1)];
    % Singular strategies at the corners (both evolving traits are at their
    % maximum or minimum values):
    elseif sum(fitgradHval<0,'all')==0 && sum(fitgradPval<0,'all')==0 && tJvalvec(end)==1 && alphavalvec(end)==finP
        tJss=[1,-10+zeros(1,maxsingstrats-1)];
        alphass=[finP,-10+zeros(1,maxsingstrats-1)];
    elseif sum(fitgradHval<0,'all')==0 && sum(fitgradPval>0,'all')==0 && tJvalvec(end)==1 && alphavalvec(1)==0
        tJss=[1,-10+zeros(1,maxsingstrats-1)];
        alphass=[0,-10+zeros(1,maxsingstrats-1)];
    elseif sum(fitgradHval>0,'all')==0 && sum(fitgradPval<0,'all')==0 && tJvalvec(1)==0 && alphavalvec(end)==finP
        tJss=[0,-10+zeros(1,maxsingstrats-1)];
        alphass=[finP,-10+zeros(1,maxsingstrats-1)];
    elseif sum(fitgradHval>0,'all')==0 && sum(fitgradPval>0,'all')==0 && tJvalvec(1)==0 && alphavalvec(1)==0
        tJss=[0,-10+zeros(1,maxsingstrats-1)];
        alphass=[0,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradHval(1,1)<0 && fitgradPval(1,1)<0 && tJvalvec(1)==0 && alphavalvec(1)==0
        tJss=[0,-10+zeros(1,maxsingstrats-1)];
        alphass=[0,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradHval(1,end)<0 && fitgradPval(1,end)>0 && tJvalvec(1)==0 && alphavalvec(end)==finP
        tJss=[0,-10+zeros(1,maxsingstrats-1)];
        alphass=[finP,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradHval(end,1)>0 && fitgradPval(end,1)<0 && tJvalvec(end)==1 && alphavalvec(1)==0
        tJss=[1,-10+zeros(1,maxsingstrats-1)];
        alphass=[0,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradHval(end,end)>0 && fitgradPval(end,end)>0 && tJvalvec(end)==1 && alphavalvec(end)==finP
        tJss=[1,-10+zeros(1,maxsingstrats-1)];
        alphass=[finP,-10+zeros(1,maxsingstrats-1)];
    elseif fitgradHval(1,1)<0 && isnan(fitgradPval(1,1)) && tJvalvec(1)==0 && alphavalvec(1)==0 % If the parasite is extinct and the host has a negative fitness gradient at (0,0)
        tJss=[0,-10+zeros(1,maxsingstrats-1)];
        alphass=[-10,-10+zeros(1,maxsingstrats-1)];
    else % No singular strategies
        tJss=-10+zeros(1,maxsingstrats);
        alphass=-10+zeros(1,maxsingstrats);
    end
end

end