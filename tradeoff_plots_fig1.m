% This code sketches different trade-off functions and loads the model 
% schematic.

%% Tolerance-reproduction trade-off

% Set up variables and parameters:
syms a(ti)
a0=5;

% Here are some example parameter values in the case of a strong,
% decelerating trade-off:
c1=0.5;
c2=2;

% Define the trade-off function:
a(ti)=a0*(1-(c1*(1-exp(-c2*ti))/(1-exp(-c2))));

% Create vectors of the values of the trade-off function:
xvec=zeros(100,1);
avec1=zeros(100,1);
for i=1:100
    xvec(i)=i/100;
    avec1(i)=a(xvec(i));
end

% Now do the same for a weak, decelerating trade-off:
c1=0.25;
c2=2;
a(ti)=a0*(1-(c1*(1-exp(-c2*ti))/(1-exp(-c2))));
avec2=zeros(100,1);
for i=1:100
    avec2(i)=a(xvec(i));
end

% Now do the same for a strong, accelerating trade-off:
c1=0.5;
c2=-1;
a(ti)=a0*(1-(c1*(1-exp(-c2*ti))/(1-exp(-c2))));
avec3=zeros(100,1);
for i=1:100
    avec3(i)=a(xvec(i));
end

% Now do the same for a weak, accelerating trade-off:
c1=0.25;
c2=-1;
a(ti)=a0*(1-(c1*(1-exp(-c2*ti))/(1-exp(-c2))));
avec4=zeros(100,1);
for i=1:100
    avec4(i)=a(xvec(i));
end

% Plot the trade-offs:
red=1/255*[215,48,39];
blue=1/255*[69,117,180];
subplot(2,3,3)
plot(xvec,avec1,'Color',red,'Linewidth',2)
hold on
plot(xvec,avec2,'Color',blue,'Linewidth',2)
hold on
plot(xvec,avec3,'Color',red,'LineStyle','--','Linewidth',2) 
hold on
plot(xvec,avec4,'Color',blue,'LineStyle','--','Linewidth',2)
xlabel('Host tolerance, $\tau_i$','interpreter','latex','fontsize',16)
ylabel('Reproduction rate, $a(\tau_i)$','interpreter','latex','fontsize',16)
ylim([2,5])
xlim([0,1])
text(0.87,4.75,'B','fontsize',24)

%%  Virulence-transmission trade-off

% Set up variables and parameters:
syms beta(alpha)
beta0=10;

% Define the trade-off function:
beta(alpha)=beta0*sqrt(alpha);

% Create vectors of the values of the trade-off function:
xvec=zeros(100,1);
betavec1=zeros(100,1);
for i=1:100
    xvec(i)=i/10;
    betavec1(i)=beta(xvec(i));
end

% Plot the trade-off:
orange=1/255*[253,174,97];
subplot(2,3,6)
plot(xvec,betavec1,'Color',orange,'Linewidth',2)
xlabel('Parasite virulence, $\alpha$','interpreter','latex','fontsize',16)
ylabel({'Transmissibility, $\beta(\alpha)$'},'interpreter','latex','fontsize',16)
ylim([0,35])
xlim([0,10])
text(8.7,32,'C','fontsize',24)

%% Include model schematic image

s=subplot(2,3,[1,2,4,5]);
A=imread('model_schematic.jpg');
B=imresize(A,2);
image(B)
set(gca,'xtick',[])
set(gca,'ytick',[])
s.Position=s.Position+[-0.01,0,0.01,0];
text(1250,50,'A','fontsize',24)
