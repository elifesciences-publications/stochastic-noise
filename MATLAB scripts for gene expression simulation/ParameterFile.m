% Defines the rate parameters used to simulate protein expression from a single allele
% takes gradient level as input to define gene ON rate and associated characteristics 
function [parameters] = ParameterFile(grade_step)
   
    %--- STEP 1 : DEFINE DEFAULT RATE CONSTANTS --------

    % Birth-death rate constants for mRNA
    Sm = 0.25;          % mRNA produced per minute from active promoter
    Dm = (0.693/15);    % mRNA decay rate - calculated from 15 min half-life measured experimentally
    
    % Birth-death rate constants for protein
    Dp = 0.693/(1*60);  % protein decay rate - 1 hour used for simulation (5 hours measured experimentally) - slowest rate constant
    expt_fano = 20;     % measured Fano factor for 2-alleles
    allele_fano = expt_fano/2;      % contribution of 1 allele
    Sp = (Dm+Dp)*(allele_fano-1);   % wildtype sens - proteins produced per minute per mRNA
    Sp = Sp*1.8;                   % microRNA site mutant sens - proteins produced per minute per mRNA
    
    % Gene ON and OFF rate constants
    Koff = 0.05;        % gene inactivation events per minute
    Kon = grade_step ;       % gene activation events per minute set along a gradient to produce a wide range of steady state protein expression
    
    % Optional - Change burst size if needed
    % mrna_burst = 4;        % Set burst size for transcription
    % Sm = mrna_burst*Koff ;  % Set Sm according to burst size 
    
    % Optional - if alleles are interacting
    ai = 1; % default - NO allelic interaction
    % fold_inhibition = 2; % Alleles inhibit each other
    % ai = 1/fold_inhibition ; % modified interaction term 
    % fold_activation = 2; % Alelles activate each other
    % ai = 1*fold_activation ; % modified interaction term  

    %--- STEP 2 : CALCULATE CHARACTERISTIC VALUES --------
    % Based on the rates in step-1, calculate steady state expression charactersitics
    Occ = Kon./(Kon+Koff);  % Promoter occupancy fraction
    Smp = Sm*Occ;           % Effective mean transcription rate
    
    % Expected Steady State mRNA and protein
    SSmRNA = Smp/Dm;            % Steady state mRNA numbers
    SSprotein = SSmRNA*Sp/Dp ;  % Steady state protein numbers
    mRNA_burst_size = Sm/Koff ; % transcription burst size
    mRNA_avg_time_between_bursts = (1/Kon) + (1/Koff);
    mRNA_burst_frequency = 1/mRNA_avg_time_between_bursts; % transcription burst frequency
    translation_burst_size_b = Sp/Dm; % translation burst size

    %--- STEP 3 : STORE VALUES IN A STRUCT CALLED parameters --------
    % rate constants
    parameters.Sm = Sm;
    parameters.Smp = Smp;
    parameters.Dm = Dm;  
    parameters.Sp = Sp ;  
    parameters.Dp = Dp;
    parameters.Kon = Kon;
    parameters.Koff = Koff;
    parameters.ai = ai;

    % expected population characateristics
    parameters.Occupancy = Occ;
    parameters.SSmRNA = SSmRNA;
    parameters.SSprotein = SSprotein;
    parameters.translation_burst_size_b = translation_burst_size_b;
    parameters.mRNA_avg_time_between_bursts = mRNA_avg_time_between_bursts ;
    parameters.mRNA_burst_frequency = mRNA_burst_frequency ;
    parameters.mRNA_burst_size = mRNA_burst_size ;