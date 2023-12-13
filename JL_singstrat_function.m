function [tolerance,virulence,host_classification,parasite_classification]=JL_singstrat_function(q,f,bJ,bA,g,tA,beta0,rJ,rA,a0,c1,c2,gamma,finP)

% This function finds the co-singular strategies and classifies their
% evolutionary and convergence stability. 

% 1 represents a co-CSS (evo stab, conv stab, strong conv stab)
% 2 is a repeller (conv unstab, strong conv unstab)
% 3 is a branching point (evo unstab, conv stab, strong conv stab)
% 4 represents an evolutionarily stable ambiguous case (evo stab, conv stab
% when the parasite evolves much faster than the host, strong conv unstab)
% 5 represents an evolutionarily unstable ambiguous case (evo unstab, conv
% stab when the parasite evolves much faster than the host, strong conv unstab)

% Set up parameters to use later:
initvec=[0.07,0.07,0.02,0.02];
orig_tmax=100;
maxsingstrats=1000;
SSres=500;

%% Section 1: Find singular strategies

% To find the singular strategies, we consider a mesh of values of tolJ
% and alpha. For each pair, we find the endemic equilibrium and then the
% two fitness gradients. We then look to see where the sign of these
% fitness gradients changes.

startH=0;
startP=0;
finH=1;
% Find the fitness gradients:
[fitgradHval,fitgradPval,tJvalvec,alphavalvec,R0counter]=JL_fitness_gradients(SSres,startH,startP,finH,finP,a0,g,q,c1,c2,rJ,rA,tA,f,beta0,bJ,bA,gamma,initvec,orig_tmax);
                                                                          
% This function determines where the fitness gradients change sign:
[markerH,markerP,~]=fitgrad_signchange_function(fitgradHval,fitgradPval);

% If no sign changes are detected:
    
% We need to consider the cases where there are no singular strategies. 
% There could still be an evolutionary attractor when one trait is at zero 
% or one:
[tJss1,alphass1]=singstrats_at_0or1(fitgradHval,fitgradPval,tJvalvec,alphavalvec,SSres,maxsingstrats,R0counter,finP);

% Back to the cases where sign changes have been found:
data=(markerH'&markerP');
data=smoothdata(smoothdata(double(data'),1,'gaussian',5),2,'gaussian',5);
data2 = data';
[~, locs1, ~, ~] = findpeaks(double(data(:))); % peaks along x
[~, locs2, ~, ~] = findpeaks(double(data2(:))); % peaks along y
data_size = size(data); % Gets matrix dimensions
[col2, row2] = ind2sub(data_size, locs2); % Converts back to 2D indices
locs2 = sub2ind(data_size, row2, col2);
ind = intersect(locs1, locs2); % Finds common peak position
[row, column] = ind2sub(data_size, ind);

% Having found the sign changes, we now 'zoom in' around them to
% determine the singular strategies more accurately:
zoomsize=10/SSres;
Hsingstratvec=-10+zeros(maxsingstrats,1);
Psingstratvec=-10+zeros(maxsingstrats,1);
for j=1:length(row)
    
    tJval=tJvalvec(row(j));
    alphaval=alphavalvec(column(j));
    
    startHtemp = max(startH,tJval - zoomsize);
    finHtemp = min(finH,tJval + zoomsize);
    startPtemp = max(startP,alphaval - zoomsize);
    finPtemp = min(finP,alphaval + zoomsize);

    % Find the fitness gradients:
    [fitgradHvalzoomed,fitgradPvalzoomed,tJvalveczoomed,alphavalveczoomed,~]=JL_fitness_gradients(SSres,startHtemp,startPtemp,finHtemp,finPtemp,a0,g,q,c1,c2,rJ,rA,tA,f,beta0,bJ,bA,gamma,initvec,orig_tmax);

    % Find where they change sign:
    [markerHzoomed,markerPzoomed,~]=fitgrad_signchange_function(fitgradHvalzoomed,fitgradPvalzoomed);
    
    data = (markerHzoomed'&markerPzoomed');
    data = smoothdata(smoothdata(double(data'),1,'gaussian',5),2,'gaussian',5);
    data2 = data';
    [~, locs1, ~, ~] = findpeaks(double(data(:))); % peaks along x
    [~, locs2, ~, ~] = findpeaks(double(data2(:))); % peaks along y
    data_size = size(data); % Gets matrix dimensions
    [col2, row2] = ind2sub(data_size, locs2); % Converts back to 2D indices
    locs2 = sub2ind(data_size, row2, col2);
    ind = intersect(locs1, locs2); % Finds common peak position
    [rowzoomed, columnzoomed] = ind2sub(data_size, ind);
    
    % If no sign changes have been detected, we look at the edges
    % of the region again:
    if isempty(rowzoomed)
        [Hsingstrats,Psingstrats]=singstrats_at_0or1(fitgradHvalzoomed,fitgradPvalzoomed,tJvalveczoomed,alphavalveczoomed,SSres,maxsingstrats,R0counter,finP);
        
        % If sign changes have been found, we record their locations and 
        % this tells us the values of the singular strategies:
        if ~any(Hsingstrats+10)
            Hsingstratvec(j)=tJval;
            Psingstratvec(j)=alphaval;
        else
            Hsingstratvec(j)=Hsingstrats(1);
            Psingstratvec(j)=Psingstrats(1);
        end
        
    else
        Hsingstratvec(j)=tJvalveczoomed(round(mean(rowzoomed)));
        Psingstratvec(j)=alphavalveczoomed(round(mean(columnzoomed)));
    end
    
end

% Create final vectors of singular strategies:
tJss=[tJss1,transpose(Hsingstratvec)];
alphass=[alphass1,transpose(Psingstratvec)];

% If alpha has evolved to the highest value allowed by the code, then we
% assume that it will evolve to increase indefinitely. 
alphass(alphass==finP)=inf;

%% Section 2: Classify Singular Strategies

% Set up vectors to be used later:
classificationH=zeros(1,maxsingstrats*2);
classificationP=zeros(1,maxsingstrats*2);

% For each singular strategy, we determine analytically whether or not it
% is evolutionarily and convergence stable:
for j=1:maxsingstrats*2
    
    if tJss(j)~=-10 && ~isnan(tJss(j))
        tJval=tJss(j);
        alphaval=alphass(j);
        
        if (alphaval==inf || alphaval==-10) && tJval==1 % If the host evolves to full tolerance then the parasite always evolves to increase virulence indefinitely.
            classH=1;
            classP=1;
        elseif (alphaval==inf || alphaval==-10) % If the parasite virulence increases indefinitely and the host does not have full tolerance, then the pathogen will go extinct.
            tJss(j)=0;
            classH=1;
            classP=NaN;
        else
            [classH,classP]=JL_classification_function(q,f,g,bJ,bA,tA,beta0,rJ,rA,a0,c1,c2,gamma,tJval,alphaval,initvec,orig_tmax);
        end
         
        classificationH(j)=classH;
        classificationP(j)=classP;
        
    end
end

% Reduce the vectors so that they only show the co-singular strategies
% which have been found:
deletevec=zeros(length(tJss),1);
for i=1:length(tJss)
    if tJss(i)==-10 && alphass(i)==-10
        deletevec(i)=i;
    elseif tJss(i)==-10
        tJss(i)=NaN;
    elseif alphass(i)==-10
        alphass(i)=NaN;
    end
end
deletevec=nonzeros(deletevec);
tJss(deletevec)=[];
alphass(deletevec)=[];

% Format vectors:
classificationH=nonzeros(classificationH);
classificationP=nonzeros(classificationP);
tolerance=tJss;
virulence=alphass;
host_classification=classificationH;
parasite_classification=classificationP;

% Format output vectors:
if length(tolerance)>1
    deletevec=zeros(length(tolerance),1);
    for i=1:length(tolerance)
        if tolerance(i)==0 || tolerance(i)==1 || virulence(i)==0
            deletevec(i)=i;
        end
    end
    deletevec=nonzeros(deletevec);
    tolerance(deletevec)=[];
    virulence(deletevec)=[];
    host_classification(deletevec)=[];
    parasite_classification(deletevec)=[];
end

end
