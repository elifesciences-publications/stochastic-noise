%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Author: Bao Wang
%       Date:  2013/11/19

%    This program performs meanshift segmentation, convex and concave anaylysis, statistics calculating.
%    This program use
%    __w_His2AvmRFP-2x-red_YanYFPattP2-2x-green_Anti-Yan-Blue_disc6_USE_QUANTIFICATION_z26.tif
%    as a single input image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function yy = ShrinkedMeanshiftSegmentation(imageFileName,channel, par1,par2,par3,par4,par5,shrink,shrinklevel,resultFileName,resultTxt,contourTxt)
      filename = imageFileName;
      filename7 = resultTxt;
      filename2 = resultFileName;
      filename8 = contourTxt;
      
      Im=imread(filename);
             Image = Im(:,:,:);  %  If you want to perform segmentaton on cropped image, use this format, 460:600,385:525
             Iorig_Red=Image(:,:,channel);% Select Red channel
             [hight width] = size(Iorig_Red(:,:)); % calculate the size of the image
             I = Iorig_Red; % Prepare I to be used for segmentation
             
             I = medfilt2(Iorig_Red); % call median filter to remove salt and pepper noise
             I = adapthisteq(I); % call adaptive histogram equalization to improve image contrast
             

                para1 = par1;
                para2 = par2;     % Set parameters 1~5, Arnau: add the meaning of the para from the notes you took
                para3 = par3;
                para4 = par4;
                para5 = par5;
  
                            [S1 labels] =  msseg(I,para3,para1,para2);   % 10, 7, 30 defalt  calls segmentation step
                             
                               S = S1; % copy S1 to S for further filtering
                 
                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mean intensity globally thresholding, remove noise
                            threshold = mean(mean(S))

                                for j=1:hight
                                    for k=1:width
                                        if (S(j,k) < threshold)
                                            S(j,k) = 0;
                                            labels(j,k) = 0;
                                        end
                                    end
                                end
                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  If Average intensity within each label< para5, then remove the label
                             maxlabel=max(max(labels));
                            label2=zeros(hight,width);
                             S2=zeros(hight,width);

                             sumI = zeros(1,maxlabel);
                             count = zeros(1,maxlabel);
                             check = zeros(1,maxlabel);


                                for j=1:hight
                                    for k=1:width
                                            temp = labels(j,k);
                                            if(temp ~= 0)
                                            sumI(1,temp)=sumI(1,temp)+double(I(j,k));
                                             count(1,temp)= count(1,temp)+1;
                                            end

                                    end
                                end
                                for i = 1:maxlabel
                                    if (( double(sumI(1,i)) / (count(1,i)) ) > para5) % Arnau, refer to your notes
                                        check(1,i) = 1;
                                    end
                                end

                                for j=1:hight
                                    for k=1:width
                                        temp = labels(j,k);
                                        if( temp ~= 0 && check(1,temp) == 1)
                                            label2(j,k)=labels(j,k);
                                                S2(j,k) = 1;
                                        end
                                    end
                                end
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            se = strel('disk',3);% Morphology filtering to smooth the edges
                            S2 = imopen(S2, se);% Now S2 contains the segmentation information


                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Concave and Convex Analysis
                            PostImage=S2; % Copy S2 to PostImage, from now on PostImage contains the segmentation result.
                            SplitImage=S2;


                            bNoNeedSplit=false;% bool for indication no need to split further
                            
                              while true
                                  if bNoNeedSplit
                                      break;
                                  else
                                      [PostImage SplitImage bNoNeedSplit]=NewSegmentBasedConvexConcave(PostImage,SplitImage,para4); % Call splitting function to split the cells
                                  end
                              end
                             
                            PostImage=SplitImage; 
                            PostImage=bwareaopen(PostImage,250);% remove particle whose size is smaller then 250 pixels
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                          
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics Analysis
                                 ImageCrop = Im(:,:,:); % read in the original RGB image
                                 ImageRed=double(Im(:,:,1));% Read in the red channel
                                 ImageGreen=double(Im(:,:,2));%Read in the green channel
                                 ImageBlue=double(Im(:,:,3));%Read in the blue channel
                                 
                                 
                                 
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% "Shrinked" segmented image 
                                if shrink == 1
                             
                                se = strel('disk',shrinklevel);        
                                PostImage = imerode(PostImage,se);
                                
                                end
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                 
                                 
                                 
                                [BinaryLabel,NumberofLabels] = bwlabel(PostImage); % calculate the total number of clusters
        
                                CountLabels = zeros(1, NumberofLabels);  % A series of parameter anouncement
                                SumIntensityR = zeros(1, NumberofLabels);
                                SumIntensityG = zeros(1, NumberofLabels);
                                SumIntensityB = zeros(1, NumberofLabels);
                                SumX = zeros(1, NumberofLabels);
                                SumY = zeros(1, NumberofLabels);
                                AvgRed = zeros(1, NumberofLabels);
                                AvgGreen = zeros(1, NumberofLabels);
                                AvgBlue = zeros(1, NumberofLabels);
                                CentX = zeros(1, NumberofLabels);
                                CentY = zeros(1, NumberofLabels);
                                DeviationsRed = zeros(1,NumberofLabels);
                                DeviationsGreen = zeros(1,NumberofLabels);
                                DeviationsBlue = zeros(1,NumberofLabels);
                                NumberofReduced = 0;
                                cellID = 1;% used to keep the valid cell number
                                
                                
                              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate centroid coordinates and averge intensities of each channel
                                for i = 1:hight
                                    for j = 1:width
                                        temp = BinaryLabel(i,j);
                                        if(temp ~= 0)
                                            CountLabels(1,temp) = CountLabels(1,temp)+1;
                                            SumIntensityR(1,temp) = SumIntensityR(1,temp)+ImageRed(i,j);
                                            SumIntensityG(1,temp) = SumIntensityG(1,temp)+ImageGreen(i,j);
                                            SumIntensityB(1,temp) = SumIntensityB(1,temp)+ImageBlue(i,j);
                                            SumX(1,temp) = SumX(1,temp)+i;
                                            SumY(1,temp) = SumY(1,temp)+j;
                                         end
                                    end
                                end
                               
                                
                               
                                slices = 10; % If Arnau uses a loop, the parameter slices should be defined in a loop
                               
                                
                                fileID8 = fopen(filename8,'w');
                                fprintf(fileID8,'%s\t%s\t%s\t%s\t\r\n','cellId','x-cor','y-cor','image_layer_measured');
                                
                                fileID7 = fopen(filename7,'w');
                                fprintf(fileID7,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\r\n','x-coord','y-coord',	'pixel_count',	'mean_green_channel','stand_dev_green',	'mean_red_channel','stand_dev_red',	'mean_blue_channel','stand_dev_blue', 'Column',	'Row',	'image_layer_measured',	  'cell_type');
                                
                                
                               for i = 1: NumberofLabels
                                   [r c] = find( BinaryLabel == i);
                                   Flag = 0;%used to let the printer not print
                                   
                                   extractSingleCell = (BinaryLabel == i);% print line by line is ok
                                   [CordinateXY,h] = contour(extractSingleCell ,[0 0],'g--');
                                   ContourX=zeros(1,1);
                                   ContourY=zeros(1,1);
                                   [rownumber columnnumber ] = size(CordinateXY);
                                   ContourX(1,1:columnnumber-1) = CordinateXY(1,2:end);
                                   ContourY(1,1:columnnumber-1) = CordinateXY(2,2:end);
                                   
                                   
                                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%remove outliers whose area size is larger than 2500
                                       if ( CountLabels(1,i) > 2500)
                                           NumberofReduced = NumberofReduced + 1;
                                           SumIntensityR(1,i) = 0;
                                           SumIntensityG(1,i) = 0;
                                           SumIntensityB(1,i) = 0;
                                           SumX(1,i) = 0;
                                           SumY(1,i) = 0;
                                           Flag = 1;
                                           for j = 1:CountLabels(1,i)
                                               BinaryLabel(r(j),c(j)) = 0;
                                           end
                                       end
                                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                                       CalculateDeviationR = zeros(1,CountLabels(1,i));
                                       CalculateDeviationG = zeros(1,CountLabels(1,i));
                                       CalculateDeviationB = zeros(1,CountLabels(1,i));
                                       
                                       for k = 1:CountLabels(1,i)
                                           CalculateDeviationR(1,k) = ImageRed(r(k),c(k));
                                           CalculateDeviationG(1,k) = ImageGreen(r(k),c(k));
                                           CalculateDeviationB(1,k) = ImageBlue(r(k),c(k));
                                       end
        
                                       DeviationsRed(1,i) = std(CalculateDeviationR);
                                       DeviationsGreen(1,i) = std(CalculateDeviationG);
                                       DeviationsBlue(1,i) = std(CalculateDeviationB);

        
                                       AvgRed(1,i) = double(SumIntensityR(1,i))/CountLabels(1,i);
                                       AvgGreen(1,i) = double(SumIntensityG(1,i))/CountLabels(1,i);
                                       AvgBlue(1,i) = double(SumIntensityB(1,i))/CountLabels(1,i);
                                       CentY(1,i) = double(SumX(1,i))/ CountLabels(1,i);
                                       CentX(1,i) = double(SumY(1,i))/ CountLabels(1,i);
                                       
                                       if (Flag ~= 1)
                                       statistics = [ CentX(1,i); CentY(1,i); CountLabels(1,i);AvgGreen(1,i);DeviationsGreen(1,i);AvgRed(1,i);DeviationsRed(1,i);AvgBlue(1,i);DeviationsBlue(1,i);-1;-1;(slices-1);-1];
                                       fprintf(fileID7,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t\r\n',statistics);
                                       statistics2 = [ cellID.*ones(1,columnnumber-1);ContourX; ContourY; (slices-1).*ones(1,columnnumber-1)];
                                       fprintf(fileID8,'%d\t%f\t%f\t%d\t\r\n',statistics2);
                                       cellID = cellID + 1;
                                       end
                               end
                               
                               fclose(fileID8);
                               fclose(fileID7);
                               
                               

			                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% saving the output segmented image 
                                f=figure('visible','off'); %f=figure; % to turn on use the second
                                imshow(ImageRed,[],'border','tight');hold on;contour(PostImage,[0 0],'b-','linewidth',1);
                                set(f,'PaperUnits','inches');
                                set(f,'PaperPosition',[1 1 19.40 19.40]); %%%% 19.4 is the size of the image devided by 100
                                set(f,'PaperPositionMode','manual');
                                print(f,'-dtiff','-r200',filename2);
                                close(f);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                
                                
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save statistics.txt statistics -ASCII ;
        
                           
%                             [r c] = size(CentX);
%                             statistics = [ CentX; CentY; CountLabels;AvgGreen;DeviationsGreen;AvgRed;DeviationsRed;AvgBlue;DeviationsBlue;-1*ones(1,c);-1*ones(1,c);(slices-1)*ones(1,c);-1*ones(1,c)];
%                             fileID = fopen(filename7,'w');
%                              fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\r\n','x-coord','y-coord',	'pixel_count',	'mean_green_channel','stand_dev_green',	'mean_red_channel','stand_dev_red',	'mean_blue_channel','stand_dev_blue', 'Column',	'Row',	'image_layer_measured',	  'cell_type');
%                             fprintf(fileID,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t\r\n',statistics);
%                             fclose(fileID);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

