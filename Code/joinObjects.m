function BW2 =  joinObjects(BW1,backgroundMask)

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

%-------- regular size check and check input arguments
[rows,cols]=size(BW1);
if ~exist('backgroundMask')     
    backgroundMask=zeros(size(BW1));        
end
%if ~exist('maxDistance')        maxDistance=rows;                       end
%BW1                                                 = imfilter(double(BW1),gaussF(3,3,1))>0.1;
%------ label input data to process objects by pairs, get major axis to determine the distances to other objects
[BW2,numInitialObjs]                                = bwlabel(BW1);

if numInitialObjs>1
    %----- only process if there are two or more objects in the region!

    statsObjects0                                       = regionprops(BW2,'Area','minorAxislength','majorAxislength','EulerNumber','BoundingBox');

    %filtG= gaussF(3,3,1);

    for countObjs=1:numInitialObjs
        %----- process of joining close objects, analyse each against its closest object
        currentObject                                   = (BW2==countObjs);
        %----- object could be empty if joined to another previously
        if ~isempty(currentObject)
            %----- generate distance map to find distance to other objects
            distanceMap                                 = bwdist(currentObject);
            %----- exclude current object from the all the objects
            BWT                                         = BW1&(BW2~=countObjs);
            %----- this superimposes all objects except current in the distance map and get the minimum value = min distance
            distToClosest                               = (min(distanceMap(BWT)));
            %----- the minimum distance determines the object to which it corresponds
            numClosestObj                               = max(BW2(distanceMap==distToClosest));
            %----- create a separate object to which comparisons will be made
            closestObject                               = (BW2==numClosestObj);
            %----- the major axis of each object will be used to determine how far away to look for objects,
            %----- distances will be dampened by log to allow small objects to be considered and prevent large ones looking too far away
            dist1                                       = 1+log(statsObjects0(numClosestObj).MajorAxisLength);
            dist2                                       = 1+log(statsObjects0(countObjs).MajorAxisLength );
            maxDistanceAccepted                         = dist1+dist2;
            %surfdat(ismember(BW2,[countObjs numClosestObj]));drawnow
            %[ statsObjects0(numClosestObj).MinorAxisLength,statsObjects0(countObjs).MinorAxisLength ]
            %dist1=min([ statsObjects0(numClosestObj).MinorAxisLength,statsObjects0(countObjs).MinorAxisLength ])*0.8;
            %dist2=max([ statsObjects0(numClosestObj).MinorAxisLength,statsObjects0(countObjs).MinorAxisLength ])*0.1;
            %[distToClosest dist1 dist2]
            %[distToClosest maxDistanceAccepted]
            %figure(1);ttt=imfilter(double(closestObject),gaussF(3,3,1));surfdat(bwmorph(ttt>0.1,'skel','inf')+2*(ttt>0.1))
            %figure(2);ttt=imfilter(double(BW2==countObjs),gaussF(3,3,1));surfdat(bwmorph(ttt>0.1,'skel','inf')+2*(ttt>0.1))

            %maxDistanceAccepted     = min(maxDistance,min(statsObjects0(numClosestObj).MinorAxisLength,statsObjects0(countObjs).MinorAxisLength));
            if (floor(distToClosest)<=ceil(maxDistanceAccepted))
                %---- join regions, first is turned into second to allow second to be analysed as well
                distanceMap2                            = (bwdist(closestObject))+distanceMap;
                %---- check that the region in between is NOT background region slightly enlarged
                regionOfConnection                      = (distanceMap2<=ceil(2+distToClosest));
                ratioOfBackground                       = mean(backgroundMask(regionOfConnection));
                if (ratioOfBackground < 0.6)
                    %----- if the region in between is NOT predominantly background then check shape

                    if ((statsObjects0(numClosestObj).EulerNumber)+(statsObjects0(countObjs).EulerNumber))>1
                        %----- accept a merger only between objects that do not contain holes
                        %----- the region of connection is reduced to tightest
                        regionOfConnection                  =  bwmorph((distanceMap2<=round(1+distToClosest+1e-5)),'fill');
                        %----- get skeletons from A, B, A+B, if the number of points of A+B grows, then do not accept


                        %reduce region analysed
                        %----- process of joining close objects, analyse each against its closest object
                        minCol                                          = max(1,  -2+ floor(statsObjects0(countObjs).BoundingBox(1)));
                        maxCol                                          = min(cols,minCol+4+statsObjects0(countObjs).BoundingBox(3));
                        minRow                                          = max(1,  -2+ floor(statsObjects0(countObjs).BoundingBox(2)));
                        maxRow                                          = min(rows,minRow+4+statsObjects0(countObjs).BoundingBox(4));
                        %----- process of joining close objects, analyse each against its closest object
                        minColClose                                     = max(1,  -2+ floor(statsObjects0(numClosestObj).BoundingBox(1)));
                        maxColClose                                     = min(cols,minColClose+4+statsObjects0(numClosestObj).BoundingBox(3));
                        minRowClose                                     = max(1,  -2+ floor(statsObjects0(numClosestObj).BoundingBox(2)));
                        maxRowClose                                     = min(rows,minRowClose+4+statsObjects0(numClosestObj).BoundingBox(4));


                        commonRegion = (regionOfConnection|closestObject|currentObject);
                        commonRegionRed = commonRegion ( min(minRow,minRowClose):max(maxRow,maxRowClose),min(minCol:minColClose):max(maxCol:maxColClose));
%%
                        skelCurrent                         = bwmorph(currentObject(minRow:maxRow,minCol:maxCol),'skel',inf);
                        [BPoints1,BPoints12,numPoints1]     = BranchPoints(skelCurrent);
                        skelJoint                           = bwmorph(commonRegionRed,'skel',inf);
                        [BPoints3,BPoints32,numPoints3]     = BranchPoints(skelJoint);
                        skelClosest                         = bwmorph(closestObject(minRowClose:maxRowClose,minColClose:maxColClose),'skel',inf);
                        [BPoints2,BPoints22,numPoints2]     = BranchPoints(skelClosest);
                        %[numPoints1+numPoints2 numPoints3]
%%                        
                        %figure(4);surfdat(currentObject+4*skelCurrent +  2*(closestObject)+3*skelClosest);drawnow
%%                          %subplot(141);surfdat(ismember(BW2,[countObjs numClosestObj]));
                          %subplot(142);surfdat(currentObject(minRow:maxRow,minCol:maxCol)+skelCurrent+BPoints1);
                          %subplot(143);surfdat(closestObject(minRowClose:maxRowClose,minColClose:maxColClose)+skelClosest+BPoints2);
                         % figure(5);surfdat(commonRegionRed+skelJoint+BPoints3)
%%
                        if (numPoints1+numPoints2)>=(numPoints3)
                            %--- if there are less branching points in the joined region, the Join
                            BW2(currentObject)              = numClosestObj;
                            %---- join region in between
                            BW2(regionOfConnection)         = numClosestObj;
                            %figure(3);ttt=imfilter(double(closestObject),gaussF(3,3,1));surfdat(bwmorph(ttt>0.1,'skel','inf')+2*(ttt>0.1))
                            newAxis                         = regionprops(bwlabel(closestObject),'MajorAxisLength');
                            statsObjects0(numClosestObj).MinorAxisLength=newAxis.MajorAxisLength;
                        end

                    end
                end
            end
        end
    end
    BW2=(BW2>0);
    BW2                  = bwmorph(BW2,'fill');
else
    BW2=BW1;
end