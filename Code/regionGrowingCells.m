function [finalCells,statsObjects3,finalCellsIm,statsObjects2,BW2, BW3,BW4,BW6,im2] = regionGrowingCells(varargin)
%function [finalCells,statsObjects3,finalCellsIm,BW6] = regionGrowingCells(dataIn)
%------- VARARGIN   :   {1} dataIn = image to be analysed it can be an image, mat file or file name
%-------                {4} im1, im1_HSV, im1_HSV_1, chrom3D_im1
%------- ARGOUT     :   {1} finalCells [rows x cols]
%-------
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
if nargin==0        %----- no data received, exit
    
        button = questdlg('No data or figure received','Select Input','Single File','Multiple Files in a Folder','Cancel','Cancel');
        if strcmp(button(1),'S')
            [sfile,pathname]=  uigetfile('*.*','Please Select Image to Process (jpg,bmp,tif,...)');
            if sfile~=  0
                dataInName=  strcat(pathname,'/',sfile);
                clear sfile pathname;
                dataIn= imread(dataInName);
            else
                disp('File not found');
                return;
            end
        end

    
  %  disp('No data or figure received');    finalCells=[];statsObjects3=[];finalCellsIm=[];BW6=[];  help  regionGrowingCells    return;
end

%% parse input
% Will check the input arguments, it can be either a name of a MAT file or
% an image, or it can be a variable or group of variables.
if nargin==1
    % If just one argument received, it can be a name or the image itself
    dataIn=varargin{1};
    if isa(dataIn,'char')
        % Image is a file name, can be .mat or .***
        if strcmp(dataIn(end-2:end),'mat')
            try
                % Image is a mat file with just an image
                % Image is a mat file with HSV data and Chrom3D datatry
                load (dataIn);
            catch
                disp('Could not read MAT file')
                finalCells=[];statsObjects3=[];finalCellsIm=[];BW6=[];  help  regionGrowingCells;    return;
            end

        else
            try
                dataIn=imread(dataIn);
            catch
                disp('Could not read Image file')
                finalCells=[];statsObjects3=[];finalCellsIm=[];BW6=[];  help  regionGrowingCells;    return;
            end
        end
    end
end

%% Regular size check and exit if incorrect dimensions
[rows,cols,levs]=size(dataIn);

if (rows==0)||(cols==0)||(levs==0)
    disp('error in input data');
    finalCells=[];statsObjects3=[]; finalCellsIm=[]; BW6=[];  BW4 =[];
    return;
end
%% Remove Shading

                    RCL                             = rows*cols*levs;
                    if (RCL)<1e6
                        subSampL=1;
                    elseif (RCL)<3e6
                        subSampL=2;
                    elseif (RCL)<12e6
                        subSampL=4;
                    else
                        subSampL=8;
                    end

%im1                     = shadingCorrection(dataIn(1:subSampL:end,1:subSampL:end,:));


                    [dataOutS,errSurfaceS]          = shadingCorrection(dataIn(1:subSampL:end,1:subSampL:end,:));
                    if subSampL==1
                        errSurfaceS2                 = errSurfaceS;
                    else
                        filtG = gaussF(3,3,1);
                        for counterL=1:levs
                            errSurfaceS2(:,:,counterL)  = imfilter(expandu(errSurfaceS(:,:,counterL),log2(subSampL)),filtG); %#ok<AGROW>
                        end
                    end
                    dataOut                         = double(dataIn)-errSurfaceS2(1:rows,1:cols,:);
                    avChannels                      = mean(mean(dataOut));
                    maxAvChannel                    = max(avChannels);
                    for counterL=1:levs
                        dataOut(:,:,counterL)       = maxAvChannel*(dataOut(:,:,counterL)/avChannels(counterL));
                    end
                    dataOut(dataOut>255)            = 255;
                    dataOut(dataOut<0)              = 0;

im1=dataOut;
clear dataOut


%im1 = dataIn;
%% To get a Whitish Background, equalise the levels of the three channels
meanLevChannels         = mean(mean(im1));
maxMeanLev              = max(meanLevChannels);
if maxMeanLev<140 ; maxMeanLev = 150; end
if maxMeanLev>160 ; maxMeanLev = 150; end

im2(:,:,1)              = maxMeanLev * im1(:,:,1)/ meanLevChannels(1);
im2(:,:,2)              = maxMeanLev * im1(:,:,2)/ meanLevChannels(2);
im2(:,:,3)              = maxMeanLev * im1(:,:,3)/ meanLevChannels(3);
im2(im2>255)            = 255;


%% Get seeds of background, blue cells and brown cells
% BackBlueBrown will transform image from RGB to HSV, get 3D chromaticity, quantise data, find HSV boundaries for segmentation
% and determine the areas of low/high saturation and do a small region growing.

[totMask,totMask0,totMask1,dataHue] = BackBlueBrown(im2);


BW1                                 = (totMask1==1);

%% Strong edges with canny's detector  NOT USED AT THE MOMENT
dataHSV                                             = rgb2hsv(im2/255);
dataSaturation                                      = dataHSV(:,:,2);
% % dataValue                                           = imfilter(dataHSV(:,:,1),gaussF(5,5,1));
% % minDataIn                  = min(dataValue(:));
% % maxDataIn                  = max(dataValue(:));
% % edgesCanny                                          = bwlabel(edge(dataValue,'canny',[],2.3));
% % surfdat(edgesCanny.*(totMask1~=3))
% % edgesCannyOverlapIn                                 = unique ((totMask1==1).*(totMask1~=3).*imdilate(edgesCanny,[1 1;1 1]));
% % edgesCannyOverlap                                   = ismember(edgesCanny,edgesCannyOverlapIn(2:end));
% % surfdat(edgesCannyOverlap+(totMask1==1)*2)
% % edgesCannyOverlap                 = bwmorph(bwmorph(edgesCannyOverlap,'spur',50),'clean');
% %     edgesCannyOverlap1                 = imclose(edgesCannyOverlap+(totMask1==1),ones(5,2));
% %     edgesCannyOverlap2                = imclose(edgesCannyOverlap+(totMask1==1),ones(2,5));
% %  edgesCannyOverlap3                = bwmorph(bwmorph(edgesCannyOverlap1|edgesCannyOverlap2,'fill'),'spur',5);
% %  surfdat(edgesCannyOverlap3)

%% Edge erosion

%BW0                         = (totMask1==1);
%BW1                         = imerode (imfill(BW0,'holes'),[0 1 0;1 1 1;0 1 0]);%ones(3));
%BW1L                        = bwlabel((BW1+BW0)==2);
%statsObjects0               = regionprops(BW1L,'Area');
%BW2                         = bwlabel(bwmorph(ismember(BW1L,find([statsObjects0.Area]>15)),'spur',5));
BW2                         = bwmorph(bwmorph(totMask1==1,'bridge'),'spur',5);

%% Join of objects


%----- JOIN those objects that are *very* close to each other Less than the minor axis of the smallest
%----- AND the gap is not background AND euler number is not increased
BW3                         = joinObjects(BW2,totMask1==2);

%% Closing of objects
%----- CLOSE objects which are open at a short edge "C" shapes
BW4                         = closeOpenObjects(BW3);

% Clean and label objects
%BW5                         = bwlabel(bwmorph(BW4,'majority'));
BW5                         = bwlabel(BW4);
statsObjects1               = regionprops(BW5,'Area');
BW6                         = bwlabel(ismember(BW5,find([statsObjects1.Area]>25)));
statsObjects2               = regionprops(BW6,'Area','minorAxislength','majorAxislength','orientation','EquivDiameter','Perimeter','Eccentricity','Centroid','BoundingBox','filledImage','filledarea','Image','EulerNumber' );

%% Splitting of objects
% Split Cells with inner holes
if isempty(statsObjects2)
    finalCells=[];    statsObjects3=-1;
else
    [NosplittedCells,splittedCells,NoHoleCells] = splitObjects(BW6,statsObjects2,rows,cols);
    [finalCells]            = bwlabel((NosplittedCells+splittedCells+NoHoleCells)>0);
    statsObjects3           = regionprops(finalCells, 'Eccentricity','Area','EulerNumber');
end

%% Final cleaning procedure, remove objects that are only partly brown

%finalCells is already labelled, just get the number of objects  (Label to identify objects uniquely)
numBrownObjects                             = max(finalCells(:));

%Area of objects has already been obtained before   (Obtain area of objects)
%statsTotMaskLab                             = regionprops(totMaskLabeled,'Area');

%%
hueTendsToBrown                             =  sum(sum((ismember(dataHue,(1:6)))&(dataSaturation<8)))/rows/cols;
%%
if (hueTendsToBrown>0.11)
    brownRange                              = ([1:4 30 31 32]);        %[1:5]
elseif (hueTendsToBrown<=0.11)&&(hueTendsToBrown>0.069)
    brownRange                              = ([1:5 30 31 32]);        %[1:5]
elseif (hueTendsToBrown<=0.069)
    brownRange                              = ([1:6 30  31 32]);   %[1:15 30 31 32]
end

%%

if numBrownObjects>0
%%
    relFactor(numBrownObjects)                  = 0;
    for counterBrObjs=1:numBrownObjects
    %generate reliability factor as % of pixels that are originally brown
        qqqq                                    = (dataHue(find(finalCells==counterBrObjs))); %#ok<FNDSB>
%        qqq2                                    = (dataSaturation(find(finalCells==counterBrObjs))); %#ok<FNDSB>
        relFactor(1,counterBrObjs)              = (sum(ismember(qqqq,brownRange)))/numel(qqqq);     %#ok<AGROW> [1:4 31 32]
        relFactor(2,counterBrObjs)              = statsObjects3(counterBrObjs).Area;
%        relFactor(3,counterBrObjs)              = (sum(qqq2>4))/numel(qqqq);     %#ok<AGROW>
    end
%%

    %discard ALL OBJECTS whose relFactor <0.2
    finalCells(ismember(finalCells,find(relFactor(1,:)<0.25)))                              = 0;
    %discard ALL OBJECTS whose area <30
    finalCells(ismember(finalCells,find([statsObjects3.Area]<30)))                         = 0;

    %discard objects whose relFactor <0.7 AND its area is less than 100 pixels
%     finalCells( ismember( finalCells,intersect((find(relFactor(1,:)<0.85)), (find([statsObjects3.Area]<100)))  ))  = 0;
%     %discard objects whose relFactor <0.55 AND its area is less than 200 pixels
%     finalCells( ismember( finalCells,intersect((find(relFactor(1,:)<0.65)), (find([statsObjects3.Area]<200)))  ))  = 0;
%     %discard objects whose relFactor <0.4 AND its area is less than 300 pixels
%     finalCells( ismember( finalCells,intersect((find(relFactor(1,:)<0.45)), (find([statsObjects3.Area]<300)))  ))  = 0;
    
%%    
    finalCells                                  = bwlabel(finalCells>0);
%%
end
%%
statsObjects3           = regionprops(finalCells, 'Eccentricity','Area','EulerNumber');




%% Obtain border cells
borderCells1                = zeros(rows,cols,3);
%borderCells1(:,:,2) = expandu(reduceu(  (zerocross(-0.5+(finalCells>0))) ))>0;
borderCells1(:,:,2)         = zerocross(-0.5+(finalCells>0)) ;
borderCells2                = repmat(borderCells1(:,:,2),[1 1 3]);
finalCellsIm                = uint8((im2).*(1-borderCells2) + 105.*(borderCells1));
figure(2);
imagesc(finalCellsIm)




%% Region growing subfunction,
% All the region growing and definition of boundary conditions are placed
% here
function [BW0,backgroundMaskDil]=regionGrowSeeds(totMask,dataHue3,hueBrowns,hueTendsToBlue,dataSaturation,dataValue)

%-------- regular size check and determination of first zoom region
[rows,cols]=size(dataHue3);

%-------- the threshold is also estimated by assuming the distributions of foreground and background
% BackgroundMask will reject everything and endoCellsMask is the seeded regions from which to grow.
backgroundMask              = (totMask==2);
endoCellsMask               = (totMask==3);

[initialSeeds,numSeeds]     = bwlabel  (bwmorph(endoCellsMask(:,:,1),'clean'),4);
fourConnectedKernel         = strel    ([0 1 0;1 1 1; 0 1 0]);
innerBoundary               =  initialSeeds-imerode  (initialSeeds,fourConnectedKernel);
outerBoundary               = -initialSeeds+imdilate (initialSeeds,fourConnectedKernel);

seedsProps                  = regionprops(initialSeeds, 'area','pixellist','pixelidxlist');
innerB_Props                = regionprops(innerBoundary,'area','pixellist','pixelidxlist'); %#ok<NASGU>
outerB_Props                = regionprops(outerBoundary,'area','pixellist','pixelidxlist');

% loop over each of the seeds

huesSeedsLocation       = max(hueBrowns);

limitNeighbourRatio     = 0.01;
% Define the adaptive thresholds according to the amount of blue received
if   (hueTendsToBlue(2)>=0.69)&(hueTendsToBlue(1)<0.035)
    limit_satValRatio       = -29;
    limit_distanceToHue     = 13;
    backgroundMaskDil           = backgroundMask;
elseif   (hueTendsToBlue(2)>=0.69)&(hueTendsToBlue(1)>=0.035)
    limit_satValRatio       = -24;
    limit_distanceToHue     = 10;
    backgroundMaskDil           = backgroundMask;
elseif  (hueTendsToBlue(2)<0.69)&(hueTendsToBlue(2)>=0.5)
    limit_satValRatio       = -24;
    limit_distanceToHue     = 10;
    backgroundMaskDil           = imdilate(bwmorph(backgroundMask,'clean'),fourConnectedKernel);
elseif  (hueTendsToBlue(2)<0.5)&(hueTendsToBlue(2)>=0.15)
    limit_satValRatio       = -20;
    limit_distanceToHue     = 8;
    backgroundMaskDil           = imdilate(bwmorph(backgroundMask,'clean'),fourConnectedKernel);
elseif   hueTendsToBlue(2)<0.15
    limit_satValRatio       = -19;
    limit_distanceToHue     = 8;
    backgroundMaskDil           = imdilate(bwmorph(backgroundMask,'clean'),fourConnectedKernel);
end


distanceToHue= [  ];
distanceToSat= [  ]; %#ok<NASGU>
distanceToVal= [  ];  %#ok<NASGU>
%distanceToWin= [  ];
ratioOfNeighboursInRange=[]; %#ok<NASGU>
ratioOfNeighboursInRange2=[]; %#ok<NASGU>
for counterSeeds=1:numSeeds
    if (sum(initialSeeds(:)==counterSeeds)>0)
        %------- define the current Seed location [hue Sat Val] and average value --
        seedLocation            = find(initialSeeds(:)==counterSeeds);
        %valsSeedsLocation       = [dataHue3(seedLocation)      dataSaturation(seedLocation)      dataValue(seedLocation)];
        %[q1,q2]                 = max(valsSeedsLocation(:,2)-valsSeedsLocation(:,3));
        %seed_Aver               = valsSeedsLocation(q2,:);

        %------- define outer Boundary
        outerBLocation           = outerB_Props(counterSeeds).PixelIdxList;
        ratioOfNeighboursInRange=1;
        k=0;
        while( ratioOfNeighboursInRange>limitNeighbourRatio)
            k=k+1;
            % First, discard all pixels that belong to the background Mask
            try
                outerBLocation(backgroundMask(outerBLocation))=[];
            catch
                %----- some pixels are outside the image
                outerBLocation(outerBLocation<1)=[];
                %ttt=1;
            end
            % Second, verify if the boundary is joining a seeded region,
            try 
                (~isempty(outerBLocation(initialSeeds(outerBLocation)>0))); %#ok<VUNUS>
                %tt=1;
            catch
                %----- some pixels are outside the image
                outerBLocation(outerBLocation>(rows*cols))=[];
            end
            if isempty (outerBLocation)
                %----- this means the region has grown until the edges of  background from all sides, end cycle.
                ratioOfNeighboursInRange=0;
            else
                if ~isempty(outerBLocation(initialSeeds(outerBLocation)>0))
                    %------- Join Regions of seeded pixels
                    overlapPixels=outerBLocation(initialSeeds(outerBLocation)>0);
                    %figure(7);                surfdat((initialSeeds==initialSeeds(overlapPixels(1)))+(initialSeeds==counterSeeds)*2)
                    %------- define outer Boundary: a dilation is necessary if an earlier object is reached
                    %newNeighbours           = outerB_Props(initialSeeds(overlapPixels(1))).PixelIdxList;
                    newNeighbours           = imdilate(initialSeeds==initialSeeds(overlapPixels(1)),fourConnectedKernel);
                    newNeighbours           = newNeighbours-(initialSeeds==initialSeeds(overlapPixels(1)));
                    %newNeighbours(newNeighbours<1)=[];
                    outerBLocation= unique( [outerBLocation;find(newNeighbours)]);

                    %change whole region
                    seedLocation = [seedLocation;seedsProps(initialSeeds(overlapPixels(1))).PixelIdxList(:,:)]; %#ok<AGROW>
                    initialSeeds(initialSeeds==initialSeeds(overlapPixels(1))) = counterSeeds;
                    outerBLocation(ismember(outerBLocation,seedLocation))=[];
                    outerBLocation(ismember(outerBLocation,overlapPixels(1)))=[];
                else
                    %------- growing process until it stops
                    outerB_Values               = ([dataHue3(outerBLocation) dataSaturation(outerBLocation)  dataValue(outerBLocation)]);
                    try
                        distanceToHue               = abs((huesSeedsLocation)-outerB_Values(:,1));
                    catch
                        %tt=1;
                    end
                    satValRatio                 =  outerB_Values(:,2)-outerB_Values(:,3);
                    distanceToSeed              = satValRatio-(distanceToHue);
                    %distanceToSat               =  abs(dataSaturation(outerBLocation)   - seed_Aver(2));
                    %distanceToVal               =  abs(dataValue(outerBLocation)        - seed_Aver(3));
                    [mindistToSeed,pixToJoin]   = max( distanceToSeed);

                    % update regions include new pixel on initialSeeds,
                    % join the pixels that is closest to the current Seed
                    initialSeeds(outerBLocation ( pixToJoin)) = counterSeeds;
                    seedsProps(counterSeeds).PixelIdxList=[seedsProps(counterSeeds).PixelIdxList;outerBLocation(pixToJoin)];
                    %ttt=1:length(distanceToHue);
                    %newMask1=zeros(rows,cols);    newMask1(outerBLocation)=distanceToHue;
                    %newMask2=zeros(rows,cols);    newMask2(outerBLocation)=satValRatio;

                    %                   figure(2);plot(distanceToWin  ,'k-d');   %surfdat(initialSeeds+3*backgroundMask);drawnow;
                    % remove pixel from outerBLocation and include in
                    % innerBLocation

                    % find new neighbours and add in outerBlocation
                    %seedsProps(counterSeeds).Area =seedsProps(counterSeeds).Area+1;
                    %----- add pixels added "outerBlocation" into the seedLocation
                    seedLocation = [seedLocation;outerBLocation(pixToJoin)]; %#ok<AGROW>
                    newNeighbours= []; %#ok<NASGU>
                    %----- find new neighbours, up, down, left and right
                    newNeighboursUp                                     = outerBLocation(pixToJoin)-1;
                    newNeighboursDown                                   = outerBLocation(pixToJoin)+1;
                    newNeighboursRight                                  = outerBLocation(pixToJoin)+rows;
                    newNeighboursLeft                                   = outerBLocation(pixToJoin)-rows;
                    %----- check that none of the neighbours are outside the image
                    boundCheck                                          = mod(outerBLocation(pixToJoin),rows);
                    newNeighboursUp    (boundCheck==1)                  =[];
                    newNeighboursDown  (boundCheck==0)                  =[];
                    newNeighboursRight (outerBLocation(pixToJoin)>(rows*(cols-1))) =[];
                    newNeighboursLeft  (outerBLocation(pixToJoin)<rows)            =[];

                    %                     if (mod(outerBLocation(pixToJoin),rows)~=1)
                    %                         newNeighbours= [newNeighbours;outerBLocation(pixToJoin)+1];
                    %                     end
                    %                     if (mod(outerBLocation(pixToJoin),rows)~=rows)
                    %                         newNeighbours= [newNeighbours;outerBLocation(pixToJoin)-1];
                    %                     end
                    %                     if (outerBLocation(pixToJoin)>rows)
                    %                         newNeighbours= [newNeighbours;outerBLocation(pixToJoin)-rows];
                    %                     end
                    %                     if (outerBLocation(pixToJoin)<(rows*(cols-1)))
                    %                         newNeighbours= [newNeighbours;outerBLocation(pixToJoin)+rows];
                    %                     end
                    %----- add new neighbours into the outerBLocation
                    outerBLocation= unique( [outerBLocation;newNeighboursRight;newNeighboursLeft;newNeighboursUp;newNeighboursDown]);
                    %outerBLocation= unique( [outerBLocation;newNeighbours]);
                    %----- remove from OuterBLocation all the pixels that belong to seedLocation

                    outerBLocation(ismember(outerBLocation,seedLocation))=[];



                    %                    ratioOfNeighboursInRange2=[ratioOfNeighboursInRange2 ;sum(((satValRatio>limit_satValRatio)&(abs(distanceToHue)<limit_distanceToHue)))/size(distanceToHue,1)];
                    ratioOfNeighboursInRange=sum(((satValRatio>limit_satValRatio)&(abs(distanceToHue)<limit_distanceToHue)))/size(distanceToHue,1);
                    %                                      if (mod(k,50)==1)
                    %                                          if ratioOfNeighboursInRange<0.01
                    %                     %                     figure(6);t=1:size(distanceToHue,1);
                    %                     %                     %plot(t,satValRatio,'b-o',t,distanceToHue,'m-d');grid on;axis tight;zoom on
                    %                     %
                    %
                    %                     %                     plot(ratioOfNeighboursInRange2,'b-d');grid on;axis tight;zoom on
                    %                     %                     %                     figure(8);imagesc(repmat((initialSeeds+backgroundMask)>0,[1 1 3]).*im1/255);drawnow;
                    %                                          figure(12);
                    %                                          %subplot(221);surfdat(newMask1);subplot(222);surfdat(newMask2);subplot(223);surfdat(newMask2-newMask1);
                    %                                          %subplot(224);
                    %                     %                     imagesc(repmat((initialSeeds+backgroundMask)>0,[1 1 3]).*im1/255);drawnow;
                    %                                            imagesc(repmat((initialSeeds)>0,[1 1 3]).*im1/255);
                    %                                            axis off;drawnow;
                    %                                            jj=jj+1
                    %                                            FRA_1(jj)=getframe;

                    %                                          %ttttt=1;
                    %                                      end

                end %finish the if - join regions or grab pixel
            end
        end% end of the loop
        %ratioOfNeighboursInRange
        %%
    end %end of the if which observes if there are any pixels in the region
    %beep;
    %figure(6);plot(ratioOfNeighboursInRange,'b-d');grid on;axis tight;zoom on
    %figure(8);%subplot(221);surfdat(newMask1);subplot(222);surfdat(newMask2);subplot(223);surfdat(newMask2-newMask1);subplot(224);
    %imagesc(repmat((initialSeeds+backgroundMask)>0,[1 1 3]).*im1/255);drawnow;
    %ttttt=1;

end
%----- clean the image from spurious branches, isolated dots and isolated holes
BW0                         = ((initialSeeds>0)&(1-backgroundMaskDil));


%BW0                         = imclose(imopen((initialSeeds>0)&(1-backgroundMaskDil),fourConnectedKernel),fourConnectedKernel);
%BW1                         = bwmorph(bwmorph(bwmorph(BW0,'hbreak'),'spur',5),'clean');
%BW1                         = bwmorph(bwmorph(bwmorph(bwmorph(initialSeeds>0,'hbreak'),'spur',5),'clean'),'fill');
%BW2                         = bwlabel(BW1);

