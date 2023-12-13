function [tolJ_end,alphaJ_end,end_pop,strain_totalH,strain_totalP,indexH_end,indexP_end,TOLJ,ALPHAJ,DISPREV,NVEC] = LL_simulation_function(t_max,a0,g,q,beta0,c1,c2,rJ,rA,tolJmin,tolJmax,tolJ_start,alphaJmin,alphaJmax,alphaJ_start,f,bJ,bA,gamma,hostmutationprob,init_pop,strain_totalH,strain_totalP,indexH_start,indexP_start,res0,nevol)

% This function runs an evolutionary simulation for a fixed number of
% timesteps.

eqtol = 1e-3;
exttol = 1e-5;

% Set up vectors to be used later:
TolJ = linspace(tolJmin,tolJmax,res0);
AlphaJ = linspace(alphaJmin,alphaJmax,res0);
TOLJ = zeros(nevol,res0);
ALPHAJ = zeros(nevol,res0);
DISPREV = zeros(nevol,1);
NVEC = zeros(nevol,1);

% Initial conditions:
tolJ_current = tolJ_start;
alphaJ_current = alphaJ_start;
indexH_current = indexH_start;
indexP_current = indexP_start;

% Each value of ievol is one evolutionary timestep.
for ievol=1:nevol
    
    % Find the ecological equilibrium:
    [~,SJ1,SA1,IJ1,IA1,~] = LL_eco_dynamics_function(t_max,a0,g,q,c1,c2,rJ,rA,tolJ_current,alphaJ_current,f,bJ,bA,beta0,gamma,eqtol,init_pop,strain_totalH,strain_totalP);
        
    % Re-format this output into a single matrix:
    SJ=zeros(strain_totalH,strain_totalP);
    SA=zeros(strain_totalH,strain_totalP);
    IJ=zeros(strain_totalH,strain_totalP);
    IA=zeros(strain_totalH,strain_totalP);
    for j=1:strain_totalH
        for k=1:strain_totalP
            SJ(j,k) = SJ1(end,j+(k-1)*strain_totalH);
            SA(j,k) = SA1(end,j+(k-1)*strain_totalH);
            IJ(j,k) = IJ1(end,j+(k-1)*strain_totalH);
            IA(j,k) = IA1(end,j+(k-1)*strain_totalH);
        end
    end
    Nhost = SA+SJ+IA+IJ;
    
    % Remove extinct classes
    Nhostrows=sum(Nhost,2);
    Nhosttotal=sum(Nhost,'all');
    
    % See if any host strains go extinct:
    extinct = (Nhostrows/Nhosttotal)<exttol;
    strain_totalH = strain_totalH-sum(extinct);
    SA(extinct,:) = [];
    SJ(extinct,:) = [];
    IA(extinct,:) = [];
    IJ(extinct,:) = [];
    indexH_current(extinct) = [];
    tolJ_current(extinct) = [];
    
    % Update N:
    Nparasite=IJ+IA;
    Nparasitecolumns=sum(Nparasite,1);
    Nparasitetotal=sum(Nparasite,'all');
    
    % See if any parasite strains go extinct:
    if Nparasitetotal<exttol
        disp("Parasite is extinct")
        strain_totalP=0;
        tolJ_end=tolJ_current;
        alphaJ_end=alphaJ_current;
        end_pop=init_pop;
        indexH_end=indexH_current;
        indexP_end=indexP_current;
        return
    end
    extinct1 = (Nparasitecolumns/Nparasitetotal)<exttol;
    strain_totalP = strain_totalP-sum(extinct1);
    if extinct1(1)==1 && ismember(0,extinct1)
        SJ(:,find(extinct1==0,1))=SJ(:,1);
        SA(:,find(extinct1==0,1))=SA(:,1);
    end
    SA(:,extinct1) = [];
    SJ(:,extinct1) = [];
    IA(:,extinct1) = [];
    IJ(:,extinct1) = [];
    Nparasite(:,extinct1) = [];
    indexP_current(extinct1) = [];
    alphaJ_current(extinct1) = [];
    
    % Update tracker
    Nhost=SJ+SA+IJ+IA;
    Nhostrows=sum(Nhost,2);
    Nparasitecolumns=sum(Nparasite,1);
    Nhosttotal=sum(Nhost,'all'); 
    Nparasitetotal=sum(Nparasite,'all'); 
    
    % Proportion of hosts of each strain
    TOLJ(ievol,indexH_current) = Nhostrows./Nhosttotal;
    % Proportion of parasites of each strain
    ALPHAJ(ievol,indexP_current) = Nparasitecolumns./Nparasitetotal;
    % Proportion of individuals who have the disease
    DISPREV(ievol) = (sum(IA,'all')+sum(IJ,'all'))/Nhosttotal;
    % Total population density:
    NVEC(ievol) = Nhosttotal;
    
    Nhostvector=zeros(strain_totalH*strain_totalP,1);
    Nparasitevector=zeros(strain_totalH*strain_totalP,1);
    for j=1:strain_totalH
        for k=1:strain_totalP
            Nhostvector(j+(k-1)*strain_totalH)=Nhost(j,k);
            Nparasitevector(j+(k-1)*strain_totalH)=Nparasite(j,k);
        end
    end
    
    % If the mutation occurs in the host resistance:
    if (rand<hostmutationprob)
        weightedprob = Nhostvector/sum(Nhostvector);
        cumsum1 = cumsum(weightedprob);
        r1 = rand*cumsum1(end);
        mutator_loc = (find(r1<cumsum1,1));
        mutator_locH=0;
        for j=1:strain_totalH
            for k=1:strain_totalP
                if mutator_loc==j+(k-1)*strain_totalH
                   mutator_locH=j;
                end
            end
        end
        mutator = indexH_current(mutator_locH);
 
        if(mutator==1) % Mutate up
            mutant = mutator+1;
        elseif(mutator==res0) % Mutate down
            mutant = mutator-1;
        else
            if(rand>0.5) % Mutate up
                mutant = mutator+1;
            else % Mutate down
                mutant = mutator-1;
            end
        end
        if(~ismember(mutant,indexH_current)) % New strain
            strain_totalH = strain_totalH+1;
            % Update vector of host trait values:
            tolJ_current_again=NaN(1,length(tolJ_current)+1);
            for i=1:length(tolJ_current)
                tolJ_current_again(1,i)=tolJ_current(i);
            end
            tolJ_current_again(1,end)=TolJ(mutant);
            tolJ_current=tolJ_current_again;
            % Update vector of resJ trait value indices:
            indexH_current_again=NaN(1,length(indexH_current)+1);
            for i=1:length(indexH_current)
                indexH_current_again(1,i)=indexH_current(i);
            end
            indexH_current_again(1,end)=mutant;
            indexH_current=indexH_current_again;
            % Add a small population with the new trait value:
            SJ_again=NaN(size(SJ,1)+1,size(SJ,2));
            for i=1:size(SJ,1)
                for j=1:size(SJ,2)
                    SJ_again(i,j)=SJ(i,j);
                end
            end
            SJ_again(end,:)=SJ(mutator_locH,:)/10;
            SJ=SJ_again;
            SA_again=NaN(size(SA,1)+1,size(SA,2));
            for i=1:size(SA,1)
                for j=1:size(SA,2)
                    SA_again(i,j)=SA(i,j);
                end
            end
            SA_again(end,:)=SA(mutator_locH,:)/10;
            SA=SA_again;
            IJ_again=NaN(size(IJ,1)+1,size(IJ,2));
            for i=1:size(IJ,1)
                for j=1:size(IJ,2)
                    IJ_again(i,j)=IJ(i,j);
                end
            end
            IJ_again(end,:)=IJ(mutator_locH,:)/10;
            IJ=IJ_again;
            IA_again=NaN(size(IA,1)+1,size(IA,2));
            for i=1:size(IA,1)
                for j=1:size(IA,2)
                    IA_again(i,j)=IA(i,j);
                end
            end
            IA_again(end,:)=IA(mutator_locH,:)/10;
            IA=IA_again;
            
        end
    
    
    % If the mutation occurs in the parasite virulence:
    else 
        weightedprob = Nparasitevector/sum(Nparasitevector);
        cumsum1 = cumsum(weightedprob);
        r1 = rand*cumsum1(end);
        mutator_loc = (find(r1<cumsum1,1));
        mutator_locP=0;
        for j=1:strain_totalH
            for k=1:strain_totalP
                if mutator_loc==j+(k-1)*strain_totalH
                   mutator_locP=k;
                end
            end
        end
        mutator = indexP_current(mutator_locP);
    
        if(mutator==1) % Mutate up
            mutant = mutator+1;
        elseif(mutator==res0) % Mutate down
            mutant = mutator-1;
        else
            if(rand>0.5) % Mutate up
                mutant = mutator+1;
            else % Mutate down
                mutant = mutator-1;
            end
        end
        if(~ismember(mutant,indexP_current)) % New strain
            strain_totalP = strain_totalP+1;
            % Update vector of resA trait values:
            alphaJ_current_again=NaN(1,length(alphaJ_current)+1);
            for i=1:length(alphaJ_current)
                alphaJ_current_again(1,i)=alphaJ_current(i);
            end
            alphaJ_current_again(1,end)=AlphaJ(mutant);
            alphaJ_current=alphaJ_current_again;
            % Update vectors of resA trait value indices:
            indexP_current_again=NaN(1,length(indexP_current)+1);
            for i=1:length(indexP_current)
                indexP_current_again(1,i)=indexP_current(i);
            end
            indexP_current_again(1,end)=mutant;
            indexP_current=indexP_current_again;
            % Add small populations with the new trait value to the
            % population:
            SJ_again=NaN(size(SJ,1),size(SJ,2)+1);
            for i=1:size(SJ,1)
                for j=1:size(SJ,2)
                    SJ_again(i,j)=SJ(i,j);
                end
            end
            SJ_again(:,end)=zeros(size(SJ,1),1);
            SJ=SJ_again;
            SA_again=NaN(size(SA,1),size(SA,2)+1);
            for i=1:size(SA,1)
                for j=1:size(SA,2)
                    SA_again(i,j)=SA(i,j);
                end
            end
            SA_again(:,end)=zeros(size(SA,1),1);
            SA=SA_again;
            IJ_again=NaN(size(IJ,1),size(IJ,2)+1);
            for i=1:size(IJ,1)
                for j=1:size(IJ,2)
                    IJ_again(i,j)=IJ(i,j);
                end
            end
            IJ_again(:,end)=IJ(:,mutator_locP)/10;
            IJ=IJ_again;
            IA_again=NaN(size(IA,1),size(IA,2)+1);
            for i=1:size(IA,1)
                for j=1:size(IA,2)
                    IA_again(i,j)=IA(i,j);
                end
            end
            IA_again(:,end)=IA(:,mutator_locP)/10;
            IA=IA_again;
            
        end
    end
    
    % Update initial conditions
    init_pop = NaN(1,4*strain_totalH*strain_totalP);
    for i=1:strain_totalH*strain_totalP
        init_pop(4*i)=IA(i);
        init_pop(4*i-1)=IJ(i);
        init_pop(4*i-2)=SA(i);
        init_pop(4*i-3)=SJ(i);
    end
end

% Create outputs:
tolJ_end=tolJ_current;
alphaJ_end=alphaJ_current;
end_pop=init_pop;
indexH_end=indexH_current;
indexP_end=indexP_current;

end