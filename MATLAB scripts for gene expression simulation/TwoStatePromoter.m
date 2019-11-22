% Function that simulates a single allele's stochastic mRNA and protein output
% Takes simulation time and rate constants as input
% outputs mRNA and protein numbers sampled at 5-minute intervals
function [Sim] = TwoStatePromoter(TotalTime,param)
    
    %--- STEP 1 : EXTRACT RATE CONSTANTS FROM PARAMETER FILE --------
    Sm = param.parameters.Sm;     
    Dm = param.parameters.Dm;  
    Sp = param.parameters.Sp;  
    Dp = param.parameters.Dp; 
    Koff = param.parameters.Koff; 
    Kon = param.parameters.Kon;

    %---- STEP 2 : DEFINE INITIAL CONDITIONS -----
    m(1) = 0 ; % Number of mRNA molecules m at time point 1 
    n(1) = 0 ; % Number of protein molecules n at time point 1 
    T(1) = 0 ;  % Initialize start time = 0 mins
    State(1) = 0 ; % Promoter state = OFF at time point 1
                   %  States: OFF = 0, ON = 1

    i = 1;  % Initialize indexing variable
    TimeElapsed =0; % simulation time elapsed in minutes

    %---- STEP 3 : SIMULATE ALLELE EXPRESSION USING GILLESPIE'S SSA----- 
    % Stochastic trajectory for output from one promoter
    % keeps track of tau,m,n and state at each step
    while TimeElapsed<TotalTime
        %---- STEP 3.1 : CALCULATE TIME STEP AFTER WHICH AN EVENT TAKES PLACE -----
        %Calculate event probabilities (deterministic rates) per unit time
        k1=Sm;        % Probability of mRNA synthesis         
        k2=m(i)*Dm;   % Probability of mRNA decay  
        k3=m(i)*Sp;   % Probability of Protein synthesis         
        k4=n(i)*Dp;   % Probability of Protein decay 
        k5=Kon;       % Probability of Inactive --> Active promoter
        k6=Koff;      % Probability of Active --> Inactive promoter
        k0=(k1+k2+k3+k4+k5+k6);% Summed Probability of ANY event occuring
        
       % Time to the next reaction stored in the vector tau, 
       % tau is sampled randomly from a decaying exponential scaled to k0
       % tau is NOT total time elapsed, just step length
        tau(i)=1/k0*log(1/rand);
       % Total time elapsed stored in vector T
        T(i+1) = sum(tau(1:i));
        TimeElapsed = T(i+1);

        %---- STEP 3.2 : DETERMINE WHICH EVENT TAKES PLACE -----
        %Now, flip a coin to determine which one of the reactions will take
        %place. Also update the 'unchanged' status of the other variables
        CoinFlip=rand;
        if CoinFlip<=k1/k0                  % event is mRNA synthesis
            if State(i) == 1                % check if promoter is Active
            m(i+1)=m(i)+1;                  % add 1 to mRNA count
            n(i+1)=n(i);
            State(i+1)=State(i);
            else                            % if promoter is inactive - no change
                m(i+1)=m(i); 
                n(i+1)=n(i);
                State(i+1)=State(i);
            end
        elseif CoinFlip<=(k1+k2)/k0         % event is mRNA decay
            m(i+1)=m(i)-1;                  % subtract 1 from mRNA count
            n(i+1)=n(i);
            State(i+1) = State(i);
        elseif CoinFlip<=(k1+k2+k3)/k0       % event is Protein synthesis
            n(i+1)=n(i)+1;                   % add 1 to protein count
            m(i+1)=m(i);
            State(i+1) = State(i);
        elseif CoinFlip<=(k1+k2+k3+k4)/k0   % event is Protein decay
            n(i+1)=n(i)-1;                  % subtract 1 from protein count
            m(i+1)=m(i);
            State(i+1) = State(i);
        elseif CoinFlip<=(k1+k2+k3+k4+k5)/k0 % event is Inactive --> Active
            m(i+1)=m(i);
            n(i+1)=n(i);
            State(i+1) = 1;                  % toggle promoter state to ON
        else                                 % event is Active --> Inactive
            m(i+1)=m(i);
            n(i+1)=n(i);
            State(i+1) = 0;                 % toggle promoter state to OFF
        end
      i = i+1; % go to next time step
    end
    
    %---- STEP 3.3 : STORE mRNA AND PROTEIN COUNTS OVER TIME AT THE END OF SIMULATION -----
    
    %% OPTIONAL : use when you want to collect all data points (memory intensive)
    % Sim.mRNA = uint16(m);
    % Sim.protein = uint16(n);
    % Sim.state = uint8(State); 
    % Sim.time = T;
    
    % Collect data every 5 mins
    % initialize storage vectors
    vecP =[]; vecM = []; vecS =[]; vecT =[];
    % sample in intervals of 5 mins
    for ii = 5:5:TotalTime %start time: interval: end time
        Index = length(T(T < ii)); 
        vecP=[vecP,n(Index)];       % Protein count at sampling time
        vecM=[vecM,m(Index)];       % mRNA count at sampling time
        vecS=[vecS,State(Index)];   % Promoter state at sampling time
        vecT=[vecT,ii];             % simulation time in minutes
    end

    % convert vectors to suitable data type and store in output struct Sim
    Sim.mRNA = uint8(vecM);
    Sim.protein = uint16(vecP);
    Sim.state = logical(vecS); 
    Sim.time = uint16(vecT);