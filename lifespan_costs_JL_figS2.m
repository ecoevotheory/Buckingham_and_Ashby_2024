% This code determines regions of parameter space where cycling,
% bistability, parasite extinction and stable monomorphism occur. 

% Define parameter values:

% scenario=1 for lifelong tolerance and scenario=2 for juvenile tolerance:
scenario=2; 
q=1;
f=1;
g=1;
rJ=0;
rA=0;
a0=5;
c2=3;
tA=0;

tol_min=0;
tol_max=1;
alpha_min=0;
t_max=100;
res0=101;
nevol=5000;
init_pop=[0.1,0.1,0.1,0.1];
finished=0;
eps=0.02;
maxsingstrats=100;

% Set up vectors of parameter values:
gammavec=[0,1];
bvec=[1,1/2,1/3,1/4,1/5,1/6,1/7,1/8,1/9,1/10];
beta0vec=[5,10,20];
c1vec=[1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
phi=32;
hostmutationprob=1/(1+phi);

% Set up vectors to use later:
classification_matrix1=NaN(length(c1vec),length(bvec));
classification_matrix2=NaN(length(c1vec),length(bvec));
classification_matrix3=NaN(length(c1vec),length(bvec));
classification_matrix4=NaN(length(c1vec),length(bvec));
classification_matrix5=NaN(length(c1vec),length(bvec));
classification_matrix6=NaN(length(c1vec),length(bvec));

% Vary recovery rate, baseline pathogen transmissibility, trade-off 
% strength and host lifespan:
for gamma_val=1:2
    disp("gamma_val = "+gamma_val)
    for beta0_val=1:3
        disp("beta0_val = "+beta0_val)
        for c1_val=1:10
            disp("c1_val = "+c1_val)
            for b_val=1:10
                disp("b_val = "+b_val)
                
                gamma=gammavec(gamma_val);
                beta0=beta0vec(beta0_val);
                c1=c1vec(c1_val);
                bJ=bvec(b_val);
                bA=bvec(b_val);
                
                if gamma==0
                    finP=10;
                elseif gamma==1
                    finP=25;
                end
                alpha_max=finP;

                %% Section 1 - Find co-singular strategies and their stability
                
                % The classifications represent:
                % 1: a co-CSS (evo stab, conv stab, strong conv stab)
                % 2: a repeller (conv unstab, strong conv unstab)
                % 3: a branching point (evo unstab, conv stab, strong conv stab)
                % 4: an evolutionarily stable ambiguous case (evo stab, conv stab
                % when parasite evolves much faster than host, strong conv unstab)
                % 5: an evolutionarily unstable ambiguous case (evo unstab, conv
                % stab when parasite evolves much faster than host, strong conv unstab)
                % 10n+1: the classification "n" reached from analytical calculations with
                % simulations revealing cycling
                % 10n+2: the classification "n" reached from analytical calculations with
                % simulations revealing branching
                % 10n+3: the classification "n" reached from analytical calculations with
                % simulations revealing an attractor
                % 10n+4: the classification "n" reached from analytical calculations with
                % simulations revealing a repeller
                % 10n+5: the classification "n" reached from analytical calculations with
                % simulations revealing parasite extinction
                
                % Find the co-singular strategies for tolerance and virulence:
                if scenario==1
                    [tol,vir,Hclass,Pclass]=LL_singstrat_function(q,f,bJ,bA,g,beta0,rJ,rA,a0,c1,c2,gamma,finP);
                elseif scenario==2
                    [tol,vir,Hclass,Pclass]=JL_singstrat_function(q,f,bJ,bA,g,tA,beta0,rJ,rA,a0,c1,c2,gamma,finP);
                end
                
                %% Section 2 - Determine whether cycling occurs around any of the co-singular strategies
                
                % Set up initial conditions:
                strain_totalH = 1;
                strain_totalP = 1;
                Tol = linspace(tol_min,tol_max,res0);
                Alpha = linspace(alpha_min,alpha_max,res0);
                
                parfor j=1:length(tol)
                    
                    if vir(j)~=inf
                        
                        tolstart=min(tol(j)+eps,tol_max);
                        alphastart=min(vir(j)+eps,alpha_max);
                        
                        % Set up Initial Conditions
                        initialH = find(Tol>=tolstart,1);
                        initialP = find(Alpha>=alphastart,1);
                        tol_start = Tol(initialH);
                        alpha_start = Alpha(initialP);
                        indexH_start = initialH;
                        indexP_start = initialP;
                        
                        % Allow both traits to evolve
                        rng(2)
                        if scenario==1
                            [~,~,~,~,~,~,~,TOL,ALPHA,~,~] = LL_simulation_function(t_max,a0,g,q,beta0,c1,c2,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
                        elseif scenario==2
                            [~,~,~,~,~,~,~,TOL,ALPHA,~,~] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
                        end
                        
                        % Find the virulence trait value at different points in time
                        virvec=NaN(nevol,1);
                        for i=1:nevol
                            ALPHAend1=ALPHA(i,:);
                            ALPHAend=[0 ALPHAend1 0];
                            [~,SSlocs]=findpeaks(ALPHAend);
                            singstrats=Alpha(SSlocs-1);
                            % If there is more than one end-point then it must be a branching
                            % point. This is denoted by "-1".
                            if length(SSlocs)>1
                                virvec(i)=-1;
                                % If there is one end-point then this is the current value of the
                                % evolving trait:
                            elseif length(SSlocs)==1
                                virvec(i)=singstrats;
                                % Otherwise, the disease has gone extinct. This is denoted by "-2".
                            elseif isempty(SSlocs)
                                virvec(i)=-2;
                            end
                        end
                        levelout=0;
                        virvec_end=virvec(end-100:end);
                        if sum(virvec_end==virvec_end(end))==length(virvec_end)
                            levelout=1;
                        end
                        
                        % Determine whether or not virulence is cycling:
                        vir_cycling=0;
                        if sum(virvec==-1)<500 && sum(virvec==-2)==0 && levelout==0
                            truncated_virvec=virvec;
                            truncated_virvec(1:1500)=[];
                            vir_start=virvec(1000);
                            if (max(truncated_virvec)-vir_start>0.2 || vir_start-min(truncated_virvec)>0.2) && sum(truncated_virvec==vir_start)~=0
                                vir_cycling=1;
                            end
                        end
                        
                        % Find the tolerance trait value at different points in time
                        tolvec=NaN(nevol,1);
                        for i=1:nevol
                            TOLend1=TOL(i,:);
                            TOLend=[0 TOLend1 0];
                            [~,SSlocs]=findpeaks(TOLend);
                            singstrats=Tol(SSlocs-1);
                            if length(SSlocs)>1
                                tolvec(i)=-1;
                                % If there is one end-point then this is the current value of the
                                % evolving trait:
                            elseif length(SSlocs)==1
                                tolvec(i)=singstrats;
                                % Otherwise, the host has gone extinct. This is denoted by "-2".
                            elseif isempty(SSlocs)
                                tolvec(i)=-2;
                            end
                        end
                        
                        % Determine whether or not tolerance is cycling:
                        tol_cycling=0;
                        if sum(tolvec==-1)<500 && sum(tolvec==-2)==0 && levelout==0
                            truncated_tolvec=tolvec;
                            truncated_tolvec(1:1500)=[];
                            tol_start=tolvec(1000);
                            if (max(truncated_tolvec)-tol_start>0.2 || tol_start-min(truncated_tolvec)>0.2) && sum(truncated_tolvec==tol_start)~=0
                                tol_cycling=1;
                            end
                        end
                        
                        if vir_cycling==1
                            Hclass(j)=Hclass(j)*10+1;
                        end
                        if tol_cycling==1
                            Pclass(j)=Pclass(j)*10+1;
                        end
                        
                    end
                end
                
                %% Section 3 - Investigate ambiguous convergence stability cases
                
                parfor j=1:length(tol)

                    if (Hclass(j)==4 || Pclass(j)==4 || Hclass(j)==5 || Pclass(j)==5) && Hclass(j)~=21 && Pclass(j)~=21 && Hclass(j)~=41 && Pclass(j)~=41 && Hclass(j)~=51 && Pclass(j)~=51
                        
                        tolstart=min(tol(j)+((tol_max-tol_min)/res0)*10,tol_max);
                        alphastart=min(vir(j)+((alpha_max-alpha_min)/res0)*10,alpha_max);
                        
                        % Set up Initial Conditions
                        initialH = find(Tol>=tolstart,1);
                        initialP = find(Alpha>=alphastart,1);
                        tol_start = Tol(initialH);
                        alpha_start = Alpha(initialP);
                        indexH_start = initialH;
                        indexP_start = initialP;
                        
                        % Allow both traits to evolve
                        rng(1)
                        if scenario==1
                            [~,~,~,~,~,~,~,TOL1,ALPHA1,~,~] = LL_simulation_function(t_max,a0,g,q,beta0,c1,c2,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
                        elseif scenario==2
                            [~,~,~,~,~,~,~,TOL1,ALPHA1,~,~] = JL_simulation_function(t_max,a0,g,q,beta0,c1,c2,tA,rJ,rA,tol_min,tol_max,tol_start,alpha_min,alpha_max,alpha_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol);
                        end
                        
                        % Find whether the virulence trait has evolved towards or away from
                        % the co-singular strategy:
                        ALPHAend1=ALPHA1(nevol,:);
                        ALPHAend=[0 ALPHAend1 0];
                        [~,SSlocs]=findpeaks(ALPHAend);
                        if length(SSlocs)>1 % suggests branching
                            Pclass(j)=10*Pclass(j)+2;
                        elseif length(SSlocs)==1
                            if abs(SSlocs-vir(j))<abs(alphastart-vir(j))
                                Pclass(j)=10*Pclass(j)+3;
                            else
                                Pclass(j)=10*Pclass(j)+4;
                            end
                        elseif isempty(SSlocs) % suggests parasite extinction
                            Pclass(j)=10*Pclass(j)+5;
                        end
                        
                        % Find whether the tolerance trait has evolved towards or away from
                        % the co-singular strategy:
                        TOLend1=TOL1(nevol,:);
                        TOLend=[0 TOLend1 0];
                        [~,SSlocs]=findpeaks(TOLend);
                        if length(SSlocs)>1 % suggests branching
                            Hclass(j)=10*Hclass(j)+2;
                        elseif length(SSlocs)==1
                            if abs(SSlocs-tol(j))<abs(tolstart-tol(j))
                                Hclass(j)=10*Hclass(j)+3;
                            else
                                Hclass(j)=10*Hclass(j)+4;
                            end
                        elseif isempty(SSlocs) % suggests parasite extinction
                            Hclass(j)=10*Hclass(j)+5;
                        end
                        
                    end
                end
                
                %% Section 4 - Output final classification
                
                final_class=NaN(length(tol),1);
                for j=1:length(tol)
                    if sum(vir==inf)==length(vir)
                        final_class(j)=[]; % extinction
                    elseif Hclass(j)==Pclass(j) && ismember(Hclass(j),[1,43])
                        final_class(j)=1; % co-CSS
                    elseif Hclass(j)==Pclass(j) && ismember(Hclass(j),[2,44,54])
                        final_class(j)=2; % repeller
                    elseif ismember(Hclass(j),[3,42,52,53]) || ismember(Pclass(j),[3,42,52,53])
                        final_class(j)=3; % branching point
                    elseif ismember(Hclass(j),[11,21,31,41,51]) || ismember(Pclass(j),[11,21,31,41,51])
                        final_class(j)=4; % cycling
                    elseif Hclass(j)==Pclass(j) && ismember(Hclass(j),[45,55])
                        final_class(j)=[]; % host extinction
                    elseif ismember(Pclass(j),[45,55])
                        final_class(j)=[]; % parasite extinction
                    elseif ismember(Hclass(j),[2,24,44,54]) || ismember(Pclass(j),[2,24,44,54])
                        final_class(j)=2; % repeller
                    else
                        disp("possible error")
                        final_class(j)=[];
                    end
                end
                
                if isempty(final_class)
                    classification=1; % EXTINCTION
                end
                if ismember(3,final_class) && ismember(4,final_class)
                    classification=234;
                elseif ismember(2,final_class) && ismember(4,final_class)
                    classification=24;
                elseif ismember(2,final_class) && ismember(3,final_class)
                    classification=34;
                elseif ismember(4,final_class) && ismember(1,final_class)
                    classification=24;
                elseif ismember(3,final_class) && ismember(1,final_class)
                    classification=34;
                elseif ismember(2,final_class)
                    classification=4;
                elseif ismember(4,final_class)
                    classification=2;
                elseif ismember(3,final_class)
                    classification=3;
                elseif ismember(1,final_class)
                    classification=5;
                end
                
                disp(classification)
                
                if gamma_val==1
                    if beta0_val==1
                        classification_matrix1(c1_val,b_val)=classification;
                    elseif beta0_val==2
                        classification_matrix3(c1_val,b_val)=classification;
                    elseif beta0_val==3
                        classification_matrix5(c1_val,b_val)=classification;
                    end
                elseif gamma_val==2
                    if beta0_val==1
                        classification_matrix2(c1_val,b_val)=classification;
                    elseif beta0_val==2
                        classification_matrix4(c1_val,b_val)=classification;
                    elseif beta0_val==3
                        classification_matrix6(c1_val,b_val)=classification;
                    end                    
                end

            end
        end
    end
end

classification_matrix1_full=classification_matrix1;
classification_matrix2_full=classification_matrix2;
classification_matrix3_full=classification_matrix3;
classification_matrix4_full=classification_matrix4;
classification_matrix5_full=classification_matrix5;
classification_matrix6_full=classification_matrix6;

classification_matrix1(classification_matrix1==234)=2;
classification_matrix1(classification_matrix1==24)=2;
classification_matrix1(classification_matrix1==34)=4;
classification_matrix2(classification_matrix2==234)=2;
classification_matrix2(classification_matrix2==24)=2;
classification_matrix2(classification_matrix2==34)=4;
classification_matrix3(classification_matrix3==234)=2;
classification_matrix3(classification_matrix3==24)=2;
classification_matrix3(classification_matrix3==34)=4;
classification_matrix4(classification_matrix4==234)=2;
classification_matrix4(classification_matrix4==24)=2;
classification_matrix4(classification_matrix4==34)=4;
classification_matrix5(classification_matrix5==234)=2;
classification_matrix5(classification_matrix5==24)=2;
classification_matrix5(classification_matrix5==34)=4;
classification_matrix6(classification_matrix6==234)=2;
classification_matrix6(classification_matrix6==24)=2;
classification_matrix6(classification_matrix6==34)=4;

%% Create plot

% Colours for plot:
blue=1/255*[69,117,181]; % co-CSS
red=1/255*[214,48,38]; % cycling
orange=1/255*[253,174,97]; % bistability
yellow=1/255*[254,224,144]; % branching
white=1/255*[255,255,255]; % extinction
mymap=[white;red;yellow;orange;blue];

% Plot the low beta0, no recovery case:
subplot(3,2,1)
imagesc(1./bvec,c1vec,classification_matrix1)
colormap(mymap)
caxis([1,5])
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
pbaspect([10,7,1])
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
pbaspect([10,7,1])
text(0.8,5.7,'B','fontsize',18)

% Plot the intermediate beta0, no recovery case:
subplot(3,2,3)
imagesc(1./bvec,c1vec,classification_matrix3)
colormap(mymap)
caxis([1,5])
ylabel('Strength of tolerance/reproduction trade-off, $c_1$','interpreter','latex')
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
pbaspect([10,7,1])
text(0.8,5.7,'C','fontsize',18)

% Plot the intermediate beta0 case with recovery:
subplot(3,2,4)
imagesc(1./bvec,c1vec,classification_matrix4)
colormap(mymap)
caxis([1,5])
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
pbaspect([10,7,1])
text(0.8,5.7,'D','fontsize',18)

% Plot the high beta0, no recovery case:
subplot(3,2,5)
imagesc(1./bvec,c1vec,classification_matrix5)
colormap(mymap)
caxis([1,5])
xlabel('Lifespan, $1/b$','interpreter','latex')
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
pbaspect([10,7,1])
text(0.8,5.7,'E','fontsize',18)

% Plot the high beta0 case with recovery:
subplot(3,2,6)
imagesc(1./bvec,c1vec,classification_matrix6)
colormap(mymap)
caxis([1,5])
set(gca,'XDir','normal')
set(gca,'YDir','normal')
set(gca,'ytick',[0.2,0.4,0.6,0.8,1]);
set(gca,'xtick',[1,4,7,10]);
pbaspect([10,7,1])
text(0.8,5.7,'F','fontsize',18)
