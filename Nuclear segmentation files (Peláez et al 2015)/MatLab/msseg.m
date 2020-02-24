
%
% Usage:
%   [S L] = msseg(I,hs,hr,M)
%    
% Inputs:
%   I  - original image in RGB or grayscale
%   hs - spatial bandwith for mean shift analysis
%   hr - range bandwidth for mean shift analysis
%   M  - minimum size of final output regions
%
% Outputs:
%   S  - segmented image
%   L  - resulting label map
%



function [S labels] = msseg(I,hs,hr,M)
  gray = 0;
  if(size(I,3)==1)
    gray = 1;
    I = repmat(I,[1 1 3]);
  end
  
  if(nargin < 4)
    hs = 10; hr = 7; M = 30;
  end
    
  [fimg labels modes regsize grad conf] = edison_wrapper(I,@RGB2Luv,...
      'SpatialBandWidth',hs,'RangeBandWidth',hr,...
      'MinimumRegionArea',M,'speedup',3);
  
  S = Luv2RGB(fimg);
  

  if(gray == 1)
    S = rgb2gray(S);
  end
