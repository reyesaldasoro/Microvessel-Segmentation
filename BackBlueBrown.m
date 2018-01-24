function [totMask,totMask0,totMask1,dataHue]=BackBlueBrown(dataRGB)
%function [totMask,totMask0,totMask1]=BackBlueBrown(dataRGB)
%
%-------- this function obtains Blue/Brown/Background masks from the HSV representation 
%-------- of and Immunohistochemistry image stained for brown/Blue (i.e. CD31 / Hematoxilyn)
%-------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------
%------             Postdoc  Sheffield University                           ----------
%------             http://carlos-reyes.staff.shef.ac.uk                    ----------
%------  11 February 2009                                 ----------------------------
%-------------------------------------------------------------------------------------
% input data:       dataHSV: an image with cells stained by immunohistochemsitry or any other colour image
%                   [sizeHue, sizeSaturation] : sizes of the histogram
% output data:      totMask0  Single Mask for which 1=initBrown, 2=initBack, 3=initBlue;
%                   totMask   smoothed version of the previous mask, un assigned pixels go to background
%                   totMask1  same as totMask but those blue pixels that touch brown pixels are
%                             considered as brown
%
%--------------------------------------------------------------------------
%
%     Copyright (C) 2012  Constantino Carlos Reyes-Aldasoro
%
%     This file is part of the RegionGrowingCells package available through www.caiman.org.uk
%
%     The RegionGrowingCells package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%
%--------------------------------------------------------------------------
%
% This m-file is part of the RegionGrowingCells package used to analyse Immunohistochemically
% stained microscopic tissues. If you find the software useful please cite the following articles:
%
%    Reyes-Aldasoro, C.C., L Williams, S Akerman, C Kanthou and G. M. Tozer, "An automatic algorithm for the
%    segmentation and morphological analysis of microvessels in immunostained histological tumour sections",
%    Journal of Microscopy, Volume 242, Issue 3, pages 262-278, June 2011. 
%
%    Reyes-Aldasoro, C.C., Griffiths,M.K., Savas, D. and Tozer, G.M., " CAIMAN: An online algorithm 
%    repository for Cancer Image Analysis" Computer Methods and Programs in Biomedicine, Volume 103, 
%    Issue 2, August 2011, Pages 97-103. 
%
%--------------------------------------------------------------------------
%
% The authors shall not be liable for any errors or responsibility for the 
% accuracy, completeness, or usefulness of any information, or method in the content, or for any 
% actions taken in reliance thereon.
%
%--------------------------------------------------------------------------



%% Parse Input
%------ no input data is received, error -------------------------
if (nargin ~=1);                                    help BackBlueBrown; 
                                                    totMask=[]; totMask1=[]; totMask0=[]; 
                                                    return;  
end

%% Pre-processing

[rows,cols,levs]                                    = size(dataRGB); %#ok<NASGU>
dataHSV                                             = rgb2hsv(dataRGB);
%----- colourHist2 will return the chromaticity histogram and will set H,S,V in different matrices
[hs_im1,chrom3D,dataHue,dataSaturation,dataValue]   = colourHist2(dataHSV,32,32,32);

%%Plot dataIn
%  figure(1)
%  if isa(dataRGB,'double') ;         imagesc(dataRGB/255);
%  else                              imagesc(dataRGB);
%  end

%%this line plots into several figures the 1D, 2D, 3D histograms and profiles, comment for long runs
%scatterChrom3D(chrom3D,69)

[sizeSaturation,sizeHue,sizeValue]          = size(chrom3D);  %#ok<NASGU>


%%  Maximum Saturation Profile P_max_S (not used at the moment ...), 99% Profile P_99_S and hue histogram

im_H_S                                      = sum(chrom3D,3);
im_H_S_cum                                  = cumsum(im_H_S);
cummulativeLevel_99                         = ceil(repmat((0.99*im_H_S_cum(end,:)),[sizeSaturation 1]));
%saturationSpread1                          = (im_H_S>0);
saturationSpread99                          = (im_H_S_cum<cummulativeLevel_99);
%P_max_S                                    = max((saturationSpread1).*(repmat((1:sizeSaturation)',[1 sizeHue])));
P_99_S                                      = max((saturationSpread99).*(repmat((1:sizeSaturation)',[1 sizeHue])));
%m_hue                                       = sum (im_H_S);
%m_hue_sat                                   = sum (im_H_S(6:end,:));
%m_hue_nonSat                                = sum (im_H_S(1:5,:));

%% Find hues biases
hueTendsToBrown                             =  sum(sum((ismember(dataHue,(1:6)))&(dataSaturation<8)))/rows/cols;
%hueTendsToBack                              =  sum(sum((ismember(dataHue,(7:14)))&(dataSaturation<8)))/rows/cols;
%hueTendsToBlue                              =  sum(sum((ismember(dataHue,(15:23)))&(dataSaturation<8)))/rows/cols;
%hueTendsToPurple                            =  sum(sum((ismember(dataHue,(24:32)))&(dataSaturation<8)))/rows/cols;

%kernels for the region growing
kernelDil1 = [0 1 0;1 1 1;0 1 0];
kernelDil3 = ones(3);
kernelDil5 = ones(5);
kernelDil7 = ones(7);


%% Define INITIAL Ranges of hues, will serve as seed to be expanded  it depends on the characteristics (blueish, redish) of the image
% Brown:        red to ocre
% Background:   yellow to cyan
disp(hueTendsToBrown);

if (hueTendsToBrown>0.55)
    brownRange                              = (2:4);
    backRange                               = (5:16);
    blueRange                               = (17:31);
    purpleRange                             = ([1 32 ]);
elseif (hueTendsToBrown<=0.55)&&(hueTendsToBrown>0.4)
    brownRange                              = (3:4);
    backRange                               = (5:16);
    blueRange                               = (17:32);
    purpleRange                             = (1:2  );
elseif (hueTendsToBrown<=0.4)&&(hueTendsToBrown>0.2)
    brownRange                              = (1:5);
    backRange                               = (6:15);
    blueRange                               = (16:30);
    purpleRange                             = (31:32 );
elseif (hueTendsToBrown<=0.2)&&(hueTendsToBrown>0.16)
    brownRange                              = ([1:4 31 32]);                %[1:5]
    backRange                               = (5:19);
    blueRange                               = (20:26);
    purpleRange                             = (27:29 );
elseif (hueTendsToBrown<=0.16)&&(hueTendsToBrown>0.11)
    brownRange                              = ([1:4 30 31 32]);        %[1:5]
    backRange                               = (5:18);               %(6:15)
    blueRange                               = (19:23);              %(16:29)
    purpleRange                             = (27:29 );             %(30:32) 
elseif (hueTendsToBrown<=0.11)&&(hueTendsToBrown>0.06)
    brownRange                              = ([1:6 30 31 32]);        %[1:5]
    backRange                               = (7:19);               %(6:15)
    blueRange                               = (20:25);              %(16:29)
    purpleRange                             = (26:29 );             %(30:32) 
elseif (hueTendsToBrown<=0.06)&&(hueTendsToBrown>0.02)
    brownRange                              = ([1:8  30 31 32]);   %[1:15 30 31 32]
    backRange                               = (9:17);              %(16:18)
    blueRange                               = (18:26);              %(19:22)
    purpleRange                             = (27:29 );             %(23:29)
elseif (hueTendsToBrown<=0.02)
    brownRange                              = ([1:15  30 31 32]);   %[1:7 32]
    backRange                               = (16:18);              %(8:19)
    blueRange                               = (19:24);              %(20:28)
    purpleRange                             = (26:29 );             %(28:31)
end
%blueRange                                  = (18:28);
    
%%

%The huesTendTo* will help define how to set the threshold for saturation
%define the Saturation Threshold as the average of P_max_S along the hues of the background
satThresholdBack                            = ceil(mean(P_99_S(backRange )))   ;
satThresholdBrown                           = 0.5*(mean(P_99_S(brownRange)))  ;
satThresholdBlue                            = 0.5*(mean(P_99_S(blueRange )))  ;
satThresholdPurple                          = 0.5*(mean(P_99_S(purpleRange))) ;

% satThresholdBrown                           = (mean(P_99_S(brownRange))) *(0.32+hueTendsToBrown/3) ;
% satThresholdBlue                            = (mean(P_99_S(blueRange ))) *(0.2+hueTendsToBlue/3) ;
% satThresholdPurple                          = (mean(P_99_S(purpleRange)))*( 0.3+hueTendsToPurple) ;


%disp([satThresholdBack satThresholdBrown satThresholdBlue]);
%totMask=[];
%return;

%%
%----------------------- Saturation  ------------------------------------------
highSatBrown                                = (dataSaturation>=satThresholdBrown);
highSatBlue                                 = (dataSaturation>=satThresholdBlue);
highSatPurple                               = (dataSaturation>=satThresholdPurple);
lowSatBack                                  = (dataSaturation<=max(2,satThresholdBack));
verylowSat                                  = (dataSaturation<=max(2,(0.7*min(P_99_S)))) ;

meanDataValue                               = mean(dataValue(:));
disp(meanDataValue)
if meanDataValue<16
    darkValue                                   = (dataValue <= (0.80*meanDataValue));
    brightValue                                 = (dataValue >  (1.05*meanDataValue));
elseif meanDataValue>19
    darkValue                                   = (dataValue <= (0.95*meanDataValue));
    brightValue                                 = (dataValue >  (0.99*meanDataValue));
else
    darkValue                                   = (dataValue <= (0.499*sizeValue));
    brightValue                                 = (dataValue >  (0.59*sizeValue));
end
medSatBrown                                 = (dataSaturation>=(max(2,satThresholdBrown/2))) ;
medSatBlue                                  = (dataSaturation>=(max(2,satThresholdBlue/2))) ;
medSatPurple                                = (dataSaturation>=(max(2,satThresholdPurple/2))) ;


%lowSatBrown                                = (dataSaturation<satThresholdBrown);
%lowSatBlue                                 = (dataSaturation<satThresholdBlue);
%highSatBack                                = (dataSaturation>=satThresholdBack);


%----------------Define seeds of regions --------------------------------------
% Above the satThreshold will define initial seeds for brown and blue
initBrown                                   = (darkValue)&  (highSatBrown) &(ismember(dataHue,brownRange));
initBlue                                    = (darkValue)&  (highSatBlue)  &(ismember(dataHue,blueRange));

initPurple                                  = (darkValue)&  (highSatPurple)&(ismember(dataHue,purpleRange));

% Below the satThreshold and restricted to hues yellow and green will define background
initBack                                    = imdilate((brightValue)&(lowSatBack)&   (ismember(dataHue,backRange)),kernelDil3);
initBack0                                   = initBack;
%%
%-----------then extend the background ---------------------------------------
%First dilate the region
boundaryRegionBack                          = imdilate(initBack,kernelDil7);
%then select very white pixels (half of the saturation threshold) exclude pixels already assigned
%Combine range and region to form the new background 
initBack                                    = initBack|((boundaryRegionBack&verylowSat)&(~(initBlue|initBrown)));

initAreas                                   = sum([initBrown(:) initBack(:) initBlue(:)])/rows/cols;

%-----------add purple region to the smallest (brown or blue) ---------------------------------------

%ALWAYS ADD PURPLE TO BROWN
%if (initAreas(1)>initAreas(3))
    %Blue is smaller, let it grow faster
%    initBlue                                = initBlue|initPurple;   
%    kernelBrown                             = kernelDil3;
%    kernelBlue                              = kernelDil5;
%else
    %Brown is smaller
    initBrown                               = initBrown|initPurple;
    kernelBrown                             = kernelDil5;
    kernelBlue                              = kernelDil1;    
%end
totMask                                     = initBrown + 2*initBack +3*initBlue;

%remAreas = 1 -sum(initAreas);
% figure(7);
% surfdat(totMask)
% figure(8); imagesc(dataRGB/255)
% figure(9); surfdat(totMask(150:200,150:200))
% figure(10); imagesc(dataRGB(150:200,150:200,:)/255)

counterGrow                                 = 1;
changeFromPrevious                          = 11;
%totMask=[];totMask0=[];
%disp(initAreas)
%return;
numGrowthCycles                             = 9;
%% Region growing two regions blue and brown, restrict to unassigned pixels and saturation levels

while ((counterGrow<numGrowthCycles)&&(changeFromPrevious>10))

    if initAreas(1)<0.2
        %-----------first extend the brown  ---------------
        %First dilate brown region
        boundaryRegionBrown                 = (conv2(double(initBrown),kernelBrown,'same'))>0;%imdilate(initBrown,kernelBrown);
        %Overlap between dilated brown region and brown range with low saturation
        brightBrown                         = (darkValue)&(medSatBrown) ;%&(ismember(dataHue,unique([(brownRange-1) (brownRange+1)])));
        %Overlap between dilated brown region and purple range with high saturation
        darkBrownPurple                     = (darkValue)&(medSatPurple);% &(ismember(dataHue,(purpleRange(1)-1:purpleRange(end)+1)));
        %Combination of regions and restriction, it is NOT part or background or brown
        combinedRegion                      = (brightBrown|darkBrownPurple)&(~(initBack|initBlue));
        initBrown                           = initBrown|(combinedRegion&boundaryRegionBrown);
    end
    if initAreas(3)<0.2

        %-----------second extend the blue, restrict it to areas not assigned before ---------------
        %First dilate blue region
        boundaryRegionBlue                  = imdilate(initBlue,kernelBlue);
        %Overlap between dilated blue region and blue range with low saturation
        brightBlue                          = (medSatBlue) &(ismember(dataHue,(blueRange(1):blueRange(end))));
        %Overlap between dilated blue region and purple range with high saturation
        darkBluePurple                      = (medSatPurple) &(ismember(dataHue,purpleRange));
        %Combination of regions and restriction, it is NOT part or background or brown
        combinedRegion                      = (brightBlue|darkBluePurple)&(~(initBack|initBrown));
        initBlue                            = initBlue|(combinedRegion&boundaryRegionBlue);
    end
    if initAreas(2)<0.2

        %First dilate back region
        boundaryRegionBack                  = (conv2(double(initBack),kernelDil1,'same'));%imdilate(initBack,kernelDil1);
        boundaryRegionLowSat                = (conv2(double(brightValue),kernelDil1,'same')>0)&(conv2(double(lowSatBack),kernelDil3,'same')>0);
        %boundaryRegionLowSat                = imdilate(brightValue,kernelDil1)&imdilate(lowSatBack,kernelDil3);
        %Combination of regions and restriction, it is NOT part or background or brown
        combinedRegion                      = (boundaryRegionBack&boundaryRegionLowSat)&(~(initBlue|initBrown));
        initBack                            = initBack|(combinedRegion);
    end
    changeFromPrevious                      = sum(totMask(:)~=(initBrown(:) + 2*initBack(:) +3*initBlue(:)));
    totMask                                 = initBrown + 2*initBack +3*initBlue;
    counterGrow                             = counterGrow+1;
end
%return;

totMask0=totMask;

%%  Dilate only the brown area, only restriction is the hue
counterGrow                                 = 1;
changeFromPrevious                          = 11;

numGrowthCycles                             = 2;

while   ((counterGrow<numGrowthCycles)&&(changeFromPrevious>10))
    %%
    %-----------third extend the brown restricted followed by blue also restricted ---------------
    %First dilate brown region
    boundaryRegionBrown                     = imdilate(initBrown,kernelDil1);
    %Overlap between dilated brown region and brown range with low saturation
    brightBrown                             = (ismember(dataHue,unique([(brownRange-1) (brownRange+1)])));
    combinedRegion                          = (brightBrown)&(~(initBack|initBlue));
    initBrown                               = initBrown|(combinedRegion&boundaryRegionBrown);
    
    changeFromPrevious                      = sum(totMask(:)~=(initBrown(:) + 2*initBack(:) +3*initBlue(:)));
    totMask                                 = initBrown + 2*initBack +3*initBlue;
    counterGrow                             = counterGrow+1;
end





%% Assignment of isolated mediumly saturated regions
% Those areas surrounded by background which have not been assigned may be left out so assign the
% bright blues and browns (not purples) but slightly eroded

initAreas                                   = sum([initBrown(:) initBack(:) initBlue(:)])/rows/cols;


% if initAreas(1)<0.2
%     if ~exist('brightBrown','var')
%         brightBrown                         = (medSatBrown) &(ismember(dataHue,brownRange(1:end)));
%     end
%     initBrown                               = initBrown|imerode(brightBrown&(~(initBack|initBlue)),kernelDil7);
% end
% 
% if initAreas(3)<0.5
%     if ~exist('brightBlue','var')
%         brightBlue                          = (medSatBlue) &(ismember(dataHue,(blueRange(1):blueRange(end))));
%     end
%     initBlue                                =  initBlue|imerode(brightBlue&(~(initBrown|initBack)),kernelDil3);
% end
%% Removal of noise (isolated pixels or pixels in pairs)

initBrownL                                  = bwlabel(initBrown,4);
initBackL                                   = bwlabel(initBack, 4);
initBlueL                                   = bwlabel(initBlue, 4);
smallBrown                                  = regionprops(initBrownL ,'Area');
smallBack                                   = regionprops(initBackL ,'Area');
smallBlue                                   = regionprops(initBlueL ,'Area');

initBrown                                   = ismember(initBrownL,find([smallBrown.Area]>=3));
initBack                                    = ismember(initBackL, find([ smallBack.Area]>=3));
initBlue                                    = ismember(initBlueL, find([ smallBlue.Area]>=3));

%initBlue                                   = bwmorph(initBlue,'clean');
%initBrown                                  = bwmorph(initBrown,'clean');
%initBack                                   = bwmorph(initBack,'clean');

%% Region growing all regions blue and brown, restrict to unassigned pixels only
counterGrow                                 = 1;
changeFromPrevious                          = 11;

%while  ((counterGrow<numGrowthCycles)&&(changeFromPrevious>10))
%     %%
%     %-----------third extend the brown restricted followed by blue also restricted ---------------
%     %First dilate brown region
%     boundaryRegionBrown                     = imdilate(initBrown,kernelDil1);
%     %Overlap between dilated brown region and brown range with low saturation
%     brightBrown                             = (ismember(dataHue,brownRange));
%     %Overlap between dilated brown region and purple range with high and low saturation
%     BrownPurple                             = (ismember(dataHue,[blueRange purpleRange]));
%     %Combination of regions and restriction, it is NOT part or background or brown
%     combinedRegion                          = (brightBrown|BrownPurple)&(~(initBack|initBlue));
%     initBrown                               = initBrown|(combinedRegion&boundaryRegionBrown);
% 
%     %First dilate blue region
%     boundaryRegionBlue                      = imdilate(initBlue,kernelDil3);
%     %Overlap between dilated blue region and blue range with low saturation
%     brightBlue                              = (ismember(dataHue,blueRange));
%     %Overlap between dilated blue region and purple range with high saturation
%     BluePurple                              =(ismember(dataHue,purpleRange));
%     %Combination of regions and restriction, it is NOT part or background or brown
%     combinedRegion                          = (brightBlue|BluePurple)&(~(initBack|initBrown));
%     initBlue                                = initBlue|(combinedRegion&boundaryRegionBlue);
% 
%     %First dilate back region
%     boundaryRegionBack                      = imdilate(initBack,kernelDil3);
%     %Combination of regions and restriction, it is NOT part or background or brown
%     combinedRegion                          = (boundaryRegionBack)&(~(initBlue|initBrown));
%     initBack                                = initBack|(combinedRegion);
% 
% 
% 
% 
% %      figure(9)
% %      surfdat(initBrown + 2*initBack +3*initBlue)
% %      drawnow;
% %      pause(0.1)
% %      initAreas = sum([initBrown(:) initBack(:) initBlue(:)])/rows/cols;
% %      remAreas = 1 -sum(initAreas);
% %      disp([initAreas remAreas])
% 
%     
%     changeFromPrevious                      = sum(totMask(:)~=(initBrown(:) + 2*initBack(:) +3*initBlue(:)));
%     totMask                                 =initBrown + 2*initBack +3*initBlue;

%end

%% add those blue and dark pixels In contact with brown, to brown
  
% % those elements not assigned to any class will be considered as background
 totMask                                    = initBrown + 2*initBack +3*initBlue;
 totMask(totMask==0)                        = 2;
 
totMask1                                    = totMask;
darkBlueNuclei                              = ((totMask==3)&(dataValue<17));


%exclude now pixels that do not touch brown cells

darkBlueNucleiL                             = bwlabel(darkBlueNuclei);
keepElem                                    = find((darkBlueNucleiL.*imdilate((totMask0==1),kernelDil1)));%[0 0 1 0 0 ;0 1 1 1 0;1 1 1 1 1 ;0 1 1 1 0; 0 0 1 0 0])));
darkNucleiNextBrown                         = ismember(darkBlueNucleiL,darkBlueNucleiL(keepElem)); %#ok<FNDSB>
totMask1( darkNucleiNextBrown)              = 1;

%% Smoothing of the final mask


initBrown                                   = (totMask1==1);
initBack                                    = (totMask1==2);
initBlue                                    = (totMask1==3);

% the combination of a closing operator and a majority smooth very nicely, BUT it also fills in the
% holes, which will be needed later on for shape analysis, therefore, keep all the holes larger than 15
% Pixels in area.
HolesBrown                                  = imfill(initBrown,'holes');
HolesBrownL                                 = bwlabel(HolesBrown-initBrown);
HolesBrownR                                 = regionprops(HolesBrownL,'area');
HolesBrownK                                 = ismember(HolesBrownL,find([HolesBrownR.Area]>10));
                             

%% Smooth with an imclose, majority and removal of small objects
initBrown                                   = imclose(bwmorph(initBrown,'spur',2),[0 1 1 0;1 1 1 1;0 1 1 0]); %ones(5));
initBrown(HolesBrownK)                      = 0;
initBrown(initBack0)                        = 0;
%initBrown                                   = bwmorph(initBrown,'spur',3);
%initBrown                                   = bwmorph(initBrown,'majority');

%%
initBrownL                                  = bwlabel(initBrown,8);
smallBrown                                  = regionprops(initBrownL ,'Area');
initBrown                                   = ismember(initBrownL,find([smallBrown.Area]>=25));
initBluePlus                                = ismember(initBrownL,find([smallBrown.Area]<25));
%%


initBlue(initBluePlus)                      = 1;
initBlue                                    = bwmorph(initBlue,'majority');
initBlue(initBrown)                         = 0;

initBack(initBrown)                         = 0;
initBack(initBlue)                          = 0;

totMask1                                    = initBrown + 2*initBack +3*initBlue;
totMask1(totMask1==0)                       = 2;
%% Final cleaning procedure, remove objects that are only partly brown

%Label to identify objects uniquelu
[totMaskLabeled,numBrownObjects]            = bwlabel(totMask1==1);

%Obtain area of objects
statsTotMaskLab                             = regionprops(totMaskLabeled,'Area');

if numBrownObjects>0
%%
    relFactor(numBrownObjects)                  = 0;
    for counterBrObjs=1:numBrownObjects
    %generate reliability factor as % of pixels that are originally brown
        qqqq                                    = (dataHue(find(totMaskLabeled==counterBrObjs))); %#ok<FNDSB>
        relFactor(counterBrObjs)                = (sum(ismember(qqqq,brownRange)))/numel(qqqq);     %#ok<AGROW>
        %valFactor(counterBrObjs)                = mean(dataSaturation(find(totMaskLabeled==counterBrObjs)));  %#ok<AGROW>
        %relFactor(counterBrObjs)                = (sum(ismember(qqqq,[1:6 31 32])))/numel(qqqq);     %#ok<AGROW>
                       
        
    end
%%


    %discard ALL OBJECTS whose relFactor <0.2
    totMask1(ismember(totMaskLabeled,find(relFactor<0.2)))                              = 3;
    
    %discard ALL OBJECTS whose area <20
    totMask1(ismember(totMaskLabeled,find([statsTotMaskLab.Area]<20)))                  = 3;

    %discard objects whose relFactor <0.75 AND its area is less than 50 pixels
    totMask1( ismember( totMaskLabeled,intersect((find(relFactor<0.5)), (find([statsTotMaskLab.Area]<50)))  ))  = 3;
    
    %discard objects whose relFactor <0.7 AND its area is less than 100 pixels
    totMask1( ismember( totMaskLabeled,intersect((find(relFactor<0.4)), (find([statsTotMaskLab.Area]<100)))  ))  = 3;
        
    %discard objects whose relFactor <0.55 AND its area is less than 200 pixels
    totMask1( ismember( totMaskLabeled,intersect((find(relFactor<0.3)), (find([statsTotMaskLab.Area]<200)))  ))  = 3;
    
end
%%

 % initBrownDil                                = imdilate(initBrown,kernelDil3);
% %initBackDil                                = imdilate(initBack,kernelDil3);
% initBlueDil                                 = imdilate(initBlue,kernelDil3);
% %%
% initBrownL                                  = bwlabel(initBrown,8);
% initBackL                                   = bwlabel(initBack, 8);
% initBlueL                                   = bwlabel(initBlue, 8);
% smallBrown                                  = regionprops(initBrownL ,'Area');
% smallBack                                   = regionprops(initBackL ,'Area');
% smallBlue                                   = regionprops(initBlueL ,'Area');
% %%
% smallBrown2                                 = ismember(initBrownL,find([smallBrown.Area]<4));
% smallBack2                                  = ismember(initBackL, find([ smallBack.Area]<4));
% smallBlue2                                  = ismember(initBlueL, find([ smallBlue.Area]<4));
% 
% %%
% initBrown(smallBrown2)                      = 0;
% initBack(smallBack2)                        = 0;
% initBlue(smallBlue2)                        = 0;
% %%
% initBrown(smallBlue2)                       = initBrownDil(smallBlue2);
% initBlue(smallBrown2)                       = initBlueDil(smallBrown2);
% %%
% 
% if (initAreas(1)>initAreas(3))
%     %Blue is smaller, let it grow first
%     initBlue(smallBack2)                    = initBlueDil(smallBack2);
%     smallBack2                              = (smallBack2)&(~initBlue);    
%     initBrown(smallBack2)                   = initBrownDil(smallBack2);
%     
% else
%     initBrown(smallBack2)                   = initBrownDil(smallBack2);
%     smallBack2                              = (smallBack2)&(~initBrown);
%     initBlue(smallBack2)                    = initBlueDil(smallBack2);
% end


%figure(6)
%surfdat(totMask1)

% % those elements not assigned to any class will be considered as background
% totMask=initBrown + 2*initBack +3*initBlue;
% totMask(totMask==0)=2;
% 
% % to smooth the mask, assign the mode of each neighbourhood to the pivot pixel
% totMask2=im2col(totMask,[3 3]);
% totMask3=mode(totMask2);
% totMask4=reshape(totMask3,rows-2,cols-2);
% totMask(2:end-1,2:end-1)=totMask4;

%figure(10) 
%surfdat(totMask)


% 
%     %%
% 
%     %%
%     %Combine range and region
%     boundaryBrown = (boundaryRegionBrown|boundaryRangeBrown);
%     %add to the existing brown ONLY if it is not part of the blue or back
%     initBrown = initBrown + boundaryBrown&(~(initBack2|initBlue));
% 
%     %-----------then extend the blue ---------------
%     %First extend regions
%     boundaryRegionBlue = imdilate(initBlue,kernelDil5);
%     %then extend ranges
%     %highSat  = highSat-1;
%     % blueRange = [min(blueRange)-1 blueRange max(blueRange)+1];
%     boundaryRangeBlue  = (highSat) &(ismember(dataHue,blueRange));
%     %Combine range and region
%     boundaryBlue = (boundaryRegionBlue|boundaryRangeBlue);
%     %add to the existing blue ONLY if it is not part of the blue or back
%     initBlue = initBlue + boundaryBlue&(~(initBack2|initBrown));
%     %%
% 
% 
% 
% 
%     % boundaryBack = imdilate(initBack,kernelDil);
%     % initBack = initBack + boundaryBack&~(initBrown|initBlue);
%     % boundaryBlue = imdilate(initBlue,kernelDil);
%     % initBlue = initBlue + boundaryBlue&~(initBack|initBrown);
%     figure(9)
%     surfdat(initBrown + 2*initBack +3*initBlue)
%     drawnow;
%     pause(0.02)
%     initAreas = sum([initBrown(:) initBack(:) initBlue(:)])/rows/cols;
%     remAreas = 1 -sum(initAreas);
%     disp([initAreas sum(initAreas)])
% 



%br=(ismember(dataHue,[1:10 30 31 32]));
%da=(ismember(dataValue,[1:21 22 23 ]));
%sa=(ismember(dataSaturation,[3 4 5:32]));

%k1=-70+80:132;
%k2=51+150:250;
%sa_da = (dataValue<(0.96098*dataSaturation+19.19939));
%figure(1);imagesc(dataRGB(k1,k2,:)/255)
%figure(2);surfdat(sa_da(k1,k2)&br(k1,k2))
