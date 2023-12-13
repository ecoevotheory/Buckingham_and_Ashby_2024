% This code determines whether different parameter sets in the lifelong 
% tolerance scenario lead to a co-CSS, full tolerance or bistability. 

% Scenario is 1 for lifelong tolerance and 2 for juvenile tolerance:
scenario=1;

% Parameter values:
q=1;
f=1;
g=1;
rJ=0;
rA=0;
a0=5;
c2=3;
tA=0;

finP=10;
tol_min=0;
tol_max=1;
alpha_min=0;
alpha_max=finP;
t_max=100;
res0=101;
nevol=5000;
init_pop=[0.1,0.1,0.1,0.1];
finished=0;
eps=0.02;
maxsingstrats=100;

bvec=[1,1/2,1/3,1/4,1/5,1/6,1/7,1/8,1/9,1/10];
beta0vec=[5,10];
c1vec=[1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
gammavec=[0,1];

% Set up vectors to use later:
classification_matrix1=NaN(length(c1vec),length(bvec));
classification_matrix2=NaN(length(c1vec),length(bvec));
classification_matrix3=NaN(length(c1vec),length(bvec));
classification_matrix4=NaN(length(c1vec),length(bvec));

% Vary the recovery rate, pathogen baseline transmissibility, strength of
% the tolerance/reproduction trade-off and host lifespan: 
for gamma_val=1:2
    for beta0_val=1:2
        for c1_val=1:10
            for b_val=1:10
                
                gamma=gammavec(gamma_val);
                beta0=beta0vec(beta0_val);
                c1=c1vec(c1_val);
                bJ=bvec(b_val);
                bA=bvec(b_val);

                % Find the co-singular strategies for tolerance and virulence:
                if scenario==1
                    [tol,vir,Hclass,Pclass]=LL_singstrat_function(q,f,bJ,bA,g,beta0,rJ,rA,a0,c1,c2,gamma,finP);
                elseif scenario==2
                    [tol,vir,Hclass,Pclass]=JL_singstrat_function(q,f,bJ,bA,g,tA,beta0,rJ,rA,a0,c1,c2,gamma,finP);
                end
                
                counter=0;
                if length(tol)>1
                    for i=1:length(tol)
                        if tol(i)==tol(1)
                            counter=counter+1;
                        end
                    end
                    if counter==length(tol)
                        tol=tol(1);
                        vir=vir(1);
                        Hclass=Hclass(1);
                        Pclass=Pclass(1);
                    end
                end
                
                if isempty(vir) || isnan(Pclass(1))
                    classification=1; % EXTINCTION
                elseif length(tol)==1 && tol(1)==1 && vir(1)==inf
                    classification=2; %FULL TOLERANCE
                elseif length(tol)==1 && (Hclass==1 || Pclass==1)
                    classification=3; %CO-CSS
                elseif length(tol)==1 && (Hclass==2 || Pclass==2)
                    classification=4; %BISTABILITY
                elseif length(tol)==1 && (Hclass==3 || Pclass==3)
                    classification=5; %BRANCHING
                elseif length(tol)==1 && Hclass==4 && Pclass==4
                    classification=3; %CO-CSS
                elseif length(tol)==1
                    classification=4; %BISTABILITY
                elseif length(tol)>1 && ~ismember(2,Hclass) && ~ismember(2,Pclass)
                    classification=3; %CO-CSS
                elseif length(tol)>1 && ~ismember(3,Hclass) && ~ismember(3,Pclass)
                    classification=4; %BISTABILITY
                elseif length(tol)>1
                    classification=4; %BISTABILITY
                end
                
                % Format output vectors:
                if gamma_val==1
                    if beta0_val==1
                        classification_matrix1(c1_val,b_val)=classification;
                    elseif beta0_val==2
                        classification_matrix3(c1_val,b_val)=classification;
                    end
                elseif gamma_val==2
                    if beta0_val==1
                        classification_matrix2(c1_val,b_val)=classification;
                    elseif beta0_val==2
                        classification_matrix4(c1_val,b_val)=classification;
                    end                    
                end
                

            end
        end
    end
end


%% Create plot

% Colours for plot:
blue=1/255*[69,117,181]; % co-CSS
orange=1/255*[253,174,97]; % bistability
yellow=1/255*[254,224,144]; % full tolerance
lightblue=1/255*[171,217,233]; % branching
white=1/255*[255,255,255]; % extinction
mymap=[white;yellow;blue;orange;lightblue];

% Plot the low beta0, no recovery case:
subplot(3,2,1)
imagesc(1./bvec,c1vec,classification_matrix1)
colormap(mymap)
caxis([1,5])
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
axis square
text(0.8,5.7,'A','fontsize',18)

% Plot the low beta0 case with recovery:
subplot(3,2,2)
imagesc(1./bvec,c1vec,classification_matrix2)
colormap(mymap)
caxis([1,5])
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
axis square
text(0.8,5.7,'B','fontsize',18)

% Plot the high beta0, no recovery case:
subplot(3,2,3)
imagesc(1./bvec,c1vec,classification_matrix3)
colormap(mymap)
caxis([1,5])
xlabel('Lifespan, $1/b$','interpreter','latex')
ylabel('Strength of tolerance/reproduction trade-off, $c_1$','interpreter','latex')
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
axis square
text(0.8,5.7,'C','fontsize',18)

% Plot the high beta0 case with recovery:
subplot(3,2,4)
imagesc(1./bvec,c1vec,classification_matrix4)
colormap(mymap)
caxis([1,5])
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
axis square
text(0.8,5.7,'D','fontsize',18)

