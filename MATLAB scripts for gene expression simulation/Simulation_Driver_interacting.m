%% Driver file which will simulate Sens expression along a gradient for cells with interacting alleles
%--- STEP 0 : SET HOME DIRECTORY FOR ATTACHED FILES --------
home_directory = pwd;
addpath(genpath(home_directory));

%--- STEP 1 : DEFINE SIMULATION SAMPLE SIZE --------
Ncells = 5000;      % number of cells simulated for each gradient level
sim_time = 5*60;    % total simulation time in minutes - 5x protein half life 

%--- STEP 2 : DEFINE A GRADIENT OF ACTIVATION --------
grade = [0.025 0.050 0.075 0.10 0.25 0.35 0.45 0.6 0.8 1.2 1.6 2 2.4 2.7 3 3.3 3.6 3.8 4 5 10];

%--- STEP 3 : SIMULATE AN EXPRESSION SPECTRUM --------
% pre-allocate memory space for faster run-time and initialize struct TF to 0
TF(numel(TFconc)).cell(Ncells).Simulation = 0  ; 
tic % calculate run time

% add required external files for parallel processor to access
parpool('AttachedFiles',{'TwoStatePromoter_interacting.m','ParameterFile.m'});

% simulation loop for each grade step
for tt=1:length(grade) 
    tt % track in Command window
    % call gene expression rate parameters for this grade 
    TF(tt).parameters = ParameterFile(grade(tt));
 
    % multiple cells are simulated independently across different processors
    parfor ii=1:Ncells % parfor opens parpool for the first time
        % store mRNA and protein numbers across sim_time per cell
        [cell(ii)] =TwoStatePromoter_interacting(sim_time, TF(tt));
    end
    
    TF(tt).cell = cell ; % All simulated cells for grade level tt are stored here
end

% save TF as a mat file
%File1 = ['test_0001.mat']
%save (File1, 'TF', '-v7.3')
delete(gcp('nocreate')) % exit parallel pool
toc
toc/60 % run time in minutes

%exit %ONLY ADD TO FILE FOR REMOTE RUNS
