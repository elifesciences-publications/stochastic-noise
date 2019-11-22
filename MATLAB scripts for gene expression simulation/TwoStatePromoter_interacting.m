% Function that simulates a single cells's stochastic mRNA and protein output for two interacting alleles
% Takes simulation time and rate constants as input
% outputs mRNA and protein numbers sampled at 5-minute intervals
function [Sim] = TwoStatePromoter_interacting(TotalTime,param)
     
    %--- STEP 1 : EXTRACT RATE CONSTANTS FROM PARAMETER FILE --------
    Sm = param.parameters.Sm;     
    Dm = param.parameters.Dm;  
    Sp = param.parameters.Sp;  
    Dp = param.parameters.Dp; 
    Koff = param.parameters.Koff; 
    Kon = param.parameters.Kon;
    ai = param.parameters.ai; % allelic interaction term

    %---- STEP 2 : DEFINE INITIAL CONDITIONS -----
    % Initial conditions for allele-1 
    m(1) = 0 ; % Number of mRNA molecules m at time point 1 
    n(1) = 0 ; % Number of protein molecules n at time point 1 
    T(1) = 0 ;  % Initialize start time = 0 mins
    State(1) = 0 ; % Promoter state = OFF at time point 1
                   %  States: OFF = 0, ON = 1

    % Initial conditions for allele-2 - denoted by suffix-b
    mb(1) = 0;
    nb(1) = 0;
    Stateb(1) = 0;

    i = 1;  % Initialize indexing variable
    TimeElapsed =0; % simulation time elapsed in minutes

    %---- STEP 3 : SIMULATE ALLELE EXPRESSION USING GILLESPIE'S SSA----- 
    % Stochastic trajectory for output from one promoter
    % keeps track of tau,m,n and state at each step
    while TimeElapsed<TotalTime
        %---- STEP 3.1 : CALCULATE TIME STEP AFTER WHICH AN EVENT TAKES PLACE -----
        % Calculate event probabilities (deterministic rates) per unit time
        k1  = Sm;
        k1b = Sm;
        k2  = m(i) *Dm;   % Probability of mRNA decay  
        k2b = mb(i)*Dm;
        k3  = m(i) *Sp;   % Probability of Protein synthesis
        k3b = mb(i)*Sp;
        k4  = n(i) *Dp;   % Probability of Protein decay 
        k4b = nb(i)*Dp;
        k5  = Kon;
        k5b = Kon;
        k6  = Koff;
        k6b = Koff;
        
        % Modify promoter activation rates based on each other's current state        
        if (State(i) == 1) % if allele-1 is ON
            k5b = Kon*ai; % Change allele-2's Kon rate
    
        elseif (Stateb(i) == 1) % if allele-2 is ON
            k5 = Kon*transvection;  % Change allele-1's Kon rate
        end
        
        k0=(k1+k2+k3+k4+k5+k6+k1b+k2b+k3b+k4b+k5b+k6b);% Summed Probability of ANY event occuring
        
    %Time to the next reaction stored in the vector tau, sampled from
    %a decaying exponential scaled to k0
    %tau is NOT total time elapsed, just step length
        tau(i)=1/k0*log(1/rand);
    % Total time elapsed stored in vector T
        T(i+1) = sum(tau(1:i));
        TimeElapsed = T(i+1);

        %---- STEP 3.2 : DETERMINE WHICH EVENT TAKES PLACE -----
        %Now, flip a coin to determine which one of the reactions will take
        %place. Also update the 'unchanged' status of the other variables
        CoinFlip=rand;
        if CoinFlip<=k1/k0    % mRNA synthesis-1
            if State(i) == 1  %if Active promoter
            m(i+1)=m(i)+1;        mb(i+1) = mb(i); 
            n(i+1)=n(i);          nb(i+1) = nb(i);
            State(i+1)=State(i);  Stateb(i+1)=Stateb(i);
            else           % if Inactive promoter
                m(i+1)=m(i);      mb(i+1) = mb(i);
                n(i+1)=n(i);      nb(i+1) = nb(i);
                State(i+1)=State(i);Stateb(i+1)=Stateb(i);
            end
            
        elseif CoinFlip<=(k1+k1b)/k0  % mRNA synthesis-2
            if Stateb(i) == 1
                m(i+1)=m(i);          mb(i+1) = mb(i)+1; 
                n(i+1)=n(i);          nb(i+1) = nb(i);
                State(i+1)=State(i);  Stateb(i+1)=Stateb(i);
            else %If inactive promoter B
                m(i+1)=m(i);         mb(i+1) = mb(i);
                n(i+1)=n(i);         nb(i+1) = nb(i);
                State(i+1)=State(i); Stateb(i+1)=Stateb(i);
            end
                
        elseif CoinFlip<=(k1+k1b+k2)/k0 % mRNA decay-1
            m(i+1)=m(i)-1;         mb(i+1) = mb(i);
            n(i+1)=n(i);           nb(i+1) = nb(i);
            State(i+1) = State(i); Stateb(i+1) = Stateb(i);
        
        elseif CoinFlip<=(k1+k1b+k2+k2b)/k0 % mRNA decay-2
            m(i+1)=m(i);           mb(i+1) = mb(i)-1;
            n(i+1)=n(i);           nb(i+1) = nb(i);
            State(i+1) = State(i); Stateb(i+1) = Stateb(i);
        
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3)/k0 % Protein synthesis-1
            n(i+1)=n(i)+1;          nb(i+1)=nb(i);
            m(i+1)=m(i);            mb(i+1)=mb(i);
            State(i+1) = State(i);  Stateb(i+1) = Stateb(i);
        
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3+k3b)/k0 % Protein synthesis-2
            n(i+1)=n(i);            nb(i+1)=nb(i)+1;
            m(i+1)=m(i);            mb(i+1)=mb(i);
            State(i+1) = State(i);  Stateb(i+1) = Stateb(i);
        
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3+k3b+k4)/k0 % Protein decay-1
            n(i+1)=n(i)-1;          nb(i+1)=nb(i); 
            m(i+1)=m(i);            mb(i+1)=mb(i);
            State(i+1) = State(i);  Stateb(i+1) = Stateb(i);
            
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3+k3b+k4+k4b)/k0 % Protein decay-2
            n(i+1)=n(i);            nb(i+1)=nb(i)-1; 
            m(i+1)=m(i);            mb(i+1)=mb(i);
            State(i+1) = State(i);  Stateb(i+1) = Stateb(i);
        
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3+k3b+k4+k4b+k5)/k0 % Inactive --> Active -1
            m(i+1)=m(i);    mb(i+1)=mb(i);
            n(i+1)=n(i);    nb(i+1)=nb(i);
            State(i+1) = 1; Stateb(i+1) = Stateb(i);
            
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3+k3b+k4+k4b+k5+k5b)/k0 % Inactive --> Active -2
            m(i+1)=m(i);            mb(i+1)=mb(i);
            n(i+1)=n(i);            nb(i+1)=nb(i);
            State(i+1) =State(i);   Stateb(i+1) = 1;
            
        elseif CoinFlip<=(k1+k1b+k2+k2b+k3+k3b+k4+k4b+k5+k5b+k6)/k0 % Active --> Inactive -1
            m(i+1)=m(i);        mb(i+1)=mb(i);
            n(i+1)=n(i);        nb(i+1)=nb(i);
            State(i+1) = 0;     Stateb(i+1) = Stateb(i);
            
        else % Active --> Inactive -2
            m(i+1)=m(i);            mb(i+1)=mb(i);
            n(i+1)=n(i);            nb(i+1)=nb(i);
            State(i+1)=State(i);    Stateb(i+1) = 0;
            
        end
    i = i+1; %next time step
end

%---- STEP 3.3 : STORE mRNA AND PROTEIN COUNTS OVER TIME AT THE END OF SIMULATION -----
%% OPTIONAL : use when you want to collect all data points (memory intensive)
%Sim.mRNA = uint16(m);       Sim.mRNAb = uint16(mb);
%Sim.protein = uint16(n);    Sim.proteinb = uint16(nb);
%Sim.state = uint8(State);   Sim.stateb = uint8(Stateb); 
%Sim.time = T;

%collect data only for every 5 mins
vecP =[]; vecM = []; vecS =[]; vecT =[];
vecPb =[]; vecMb = []; vecSb =[];
for ii = 5:5:TotalTime %start time: interval: end time
    Index = length(T(T < ii));
    vecP=[vecP,n(Index)];
    vecM=[vecM,m(Index)];
    vecS=[vecS,State(Index)];
    vecT=[vecT,ii];
    vecPb=[vecPb,nb(Index)];
    vecMb=[vecMb,mb(Index)];
    vecSb=[vecSb,Stateb(Index)];
end
%New timepoint dataset
Sim.mRNA = uint8(vecM);         Sim.mRNAb = uint8(vecMb);
Sim.protein = uint16(vecP);     Sim.proteinb = uint16(vecPb);
Sim.state = logical(vecS);      Sim.stateb = logical(vecSb);
Sim.time = uint16(vecT);
