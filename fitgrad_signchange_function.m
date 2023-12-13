function [markerH,markerP,numberofsingstrats]=fitgrad_signchange_function(fitgradHval,fitgradPval)

% This function finds where the fitness gradients change sign.

n = length(fitgradHval);
markerH = zeros(n,n,2);
markerP = zeros(n,n,2);
for k=1:4
    
    % First, we consider sign changes in fitgradHval and then we repeat for
    % fitgradPval:
    if(k<3)
        input = sign(fitgradHval);
    else
        input = sign(fitgradPval);
    end
    
    % We surround the matrix of fitness gradient signs with either 1's or
    % -1's
    INPUT = zeros(n+2);
    if(rem(k,2)==0)
        input(isnan(input))=-1;
        INPUT(1,:)=-1;
        INPUT(n+2,:)=-1;
        INPUT(:,1)=-1;
        INPUT(:,n+2)=-1;
    else
        input(isnan(input))=1;
        INPUT(1,:)=1;
        INPUT(n+2,:)=1;
        INPUT(:,1)=1;
        INPUT(:,n+2)=1;
    end
    
    % We determine whether the signs are the same above, below, to the left
    % of and to the right of each entry:
    INPUT(2:(n+1),2:(n+1)) = input;
    UP = [INPUT(2:(n+2),:);zeros(1,n+2)];
    DOWN = [zeros(1,n+2);INPUT(1:(n+1),:)];
    LEFT = [INPUT(:,2:(n+2)),zeros(n+2,1)];
    RIGHT = [zeros(n+2,1),INPUT(:,1:(n+1))];
    A = INPUT==UP;
    B = INPUT==DOWN;
    C = INPUT==LEFT;
    D = INPUT==RIGHT;
    
    % If they are not all the same then we have detected a sign change:
    OUTPUT = (A&B&C&D);
    output = 1-OUTPUT(2:(end-1),2:(end-1));
    
    % The places where these sign changes occur give our output:
    if(k<3)
        markerH(:,:,k) = output;
    else
        markerP(:,:,k-2) = output;
    end
end

% This method also wrongly detects sign changes at the edge of the fitness
% gradient matrix. For this reason, we look for sign changes both when the
% edge is assumed to be positive and when it is assumed to be negative, and
% require that both give a sign change:
markerH = double(markerH(:,:,1)&markerH(:,:,2));
markerP = double(markerP(:,:,1)&markerP(:,:,2));

% We do not want changes from positive or negative to NaN (when the host is
% not viable) to count as sign changes, so we remove them:
markerH(isnan(fitgradHval)|isnan(fitgradPval)) = 0;
markerP(isnan(fitgradHval)|isnan(fitgradPval)) = 0;

% We determine where both fitgradJval and fitgradAval change sign:
list = find(markerP&markerH);
list(isnan(fitgradHval(list)) | isnan(fitgradPval(list)))=[];
[singstratpointerH,~] = ind2sub(size(markerP),list);

% We then know how many potential singular strategies have been found:
numberofsingstrats=length(singstratpointerH);

end
