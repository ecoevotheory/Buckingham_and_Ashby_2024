function [SJval,SAval,IJval,IAval,problem_marker]=endemic_equilibrium_function(tolJ,tA,rJ,rA,g,a,q,beta0,bJ,bA,f,alpha,gamma,initvec,orig_tmax)

% This function determines the endemic equilibrium of the ecological
% system by running the system for a long time. 

% Apply trade-offs:
beta=beta0*sqrt(alpha);
betaJ=beta*(1-rJ);
betaA=beta*(1-rA);

% Set up parameters to use later:
tol = 1e-3;
problem_marker=0;

% Run the ODE solver:
[t,y] = ode_function(orig_tmax,a,q,f,g,bJ,bA,gamma,betaJ,betaA,alpha,tolJ,tA,initvec);

% Check that the system has actually reached equilibrium:
t1=find(t>t(end)*0.9,1);

% If the system has not reached equilibrium, then we need to extend the
% time:
if(any(range(y(t1:end,:))>tol))
    [t,y] = ode_function(orig_tmax*9,a,q,f,g,bJ,bA,gamma,betaJ,betaA,alpha,tolJ,tA,y(end,:));
    
    t=t+orig_tmax;
    
    % Check again and extend again...
    t2=find(t>t(end)*0.9,1);
    if(any(range(y(t2:end,:))>tol)) 
        [t,y] = ode_function(orig_tmax*90,a,q,f,g,bJ,bA,gamma,betaJ,betaA,alpha,tolJ,tA,y(end,:));
        t=t+orig_tmax*10;
        t3=find(t>t(end)*0.9,1);
        
        % Check again and extend again...
        if(any(range(y(t3:end,:))>tol))
            [t,y] = ode_function(orig_tmax*900,a,q,f,g,bJ,bA,gamma,betaJ,betaA,alpha,tolJ,tA,y(end,:));
            t=t+orig_tmax*100;
            t4=find(t>t(end)*0.9,1);
            
            % Display an error if the system has still not reached
            % equilibrium:
            if(any(range(y(t4:end,:))>50*tol))
                problem_marker=1;
                disp("In endemic_equilibrium_function equilibrium is not reached")
            end
        end
    end
end

% The output of this function is the equilibrium state of the system:
SJval=y(end,1);
SAval=y(end,2);
IJval=y(end,3);
IAval=y(end,4);

end