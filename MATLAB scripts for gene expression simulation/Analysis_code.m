% Data analysis code to calculate Fano factor from allele simulations

%---- STEP 1 : LOAD DATA FILE ----------------
% File1 = ['/Documents/MATLAB/Sens_two_color_expression_sim_0001.mat'] % file path
% TF= load(File1);  % data stored in struct TF

%---- STEP 2 : EXTRACT SINGLE CELL ALLELE-WISE OUTPUT AT DESIRED TIME-POINT ----------------
 Time = 5*60 ;   % Define sampling time, kept at 5 hours, endpoint of simulation
 Index = Time/5; % Index number corresponding to sampling time (sample times are 5 mins apart) 

 % define 2 empty vectors to store allele-wise expression data
 vec1=[]; 
 vec2=[];
 % extract total number of cells simulated from data file
 Ncells = numel(TF(1).cell);

 % define number of bins to divide expression data
 totbins = 30 ; 

 for tt = 1:numel(TF) % for each grade of expression
    for jj=1:Ncells   % for each cell in this grade
        % append single cell protein count for allele-1 at desired time
        vec1=[vec1,TF(tt).cell(jj).allele(1).protein(Index)];
        % append single cell protein count for allele-2 at desired time
        vec2=[vec2,TF(tt).cell(jj).allele(2).protein(Index)];
    end
end

 % convert to double precision for later computations
 % equivalent to 2-color experimental data
 vec1 = double(vec1);
 vec2 = double(vec2);
 
%---- STEP 3 : BIN CELLS BY TOTAL PROTEIN LEVEL  -------
 % normalize by mean for ease of binning across datasets
 scaled_vec1 = vec1/mean(vec1); 
 scaled_vec2 = vec2/mean(vec2);
 avg_allele_exp=0.5*(scaled_vec1+scaled_vec2); % keep mean total normalized to 1
 
 % calculate bin boundaries for normalized expression
 [n,bins]=hist(avg_allele_exp,totbins);
 % sort cells into bins according to expression level
 for ii=1:length(bins)-2 % -2 to leave out last bin if low sample size, otherwise -1 
    bin(ii).ind=find(bins(ii)<avg_allele_exp & avg_allele_exp<bins(ii+1));
 end

%---- STEP 4 : CALCULATE FANO FACTOR FOR EACH BIN  -------
  for ii= 1:length(bins) -2 % for each bin except the last one
    % collect allele-wise expression data for cells in this bin
    % each cell is defined by it's index or position in the vector
    clear x y ci % delete previous versions of iterables

    x = vec1(bin(ii).ind); 
    y = vec2(bin(ii).ind);

    bin_fano(ii) = FanoFactor(x,y); % Fano factor for this bin
    ci = bootci(1000, @FanoFactor,x,y); % 95% confidence intervals for estimated Fano factor
    lower_ci(ii) = bin_fano(ii) - ci(1); % yneg
    upper_ci(ii) = ci(2) - bin_fano(ii); % ypos

    bin_protein(ii) = mean(x+y); % Average total protein expression in binned cells
  end

  set(0,'DefaultFigureWindowStyle','docked') % toggle to dock
  
  figure(1) % Fano factor plot
  set(gca,'fontsize', 18);
  errorbar(bin_protein(1:end-5), bin_fano(1:end-5), lower_ci(1:end-5), upper_ci(1:end-5), 'LineWidth',2)
  title('Fano Factor with 95% CI vs Mean total protein')
  xlim([0 1500])
  ylim([0 150])
  xlabel('Sens protein (molecules)')
  ylabel('Fano factor (molecules)')
  hold on

  figure(2) % Scatter plot
  set(gca,'fontsize', 18);
  plot(vec1,vec2,'.')
  title('Protein output correlation between alleles')
  xlabel('Allele 1 protein (molecules)');
  ylabel('Allele 2 protein (molecules)'); 
  hold on;