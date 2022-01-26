function [NosplittedCells,splittedCells,NoHoleCells]= splitObjects(largeEndoCells,statsObjects,rows,cols)
%function [NosplittedCells,splittedCells,NoHoleCells]= splitObjects(largeEndoCells,statsObjects,rows,cols)
%
%-------- this function processes a series of objects (endothelial cells) and splits into two or more objects
%-------- those objects with two or more inner holes (lumen) which are assumed to be joint vessels
%-------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------
%------             Postdoc  Sheffield University                           ----------
%------             http://tumour-microcirculation.group.shef.ac.uk         ----------
%------  5 June    2008                                                     ----------
%-------------------------------------------------------------------------------------
% input data:       largeEndoCells  : the image with the labelled Objects
%                   statsObjects    : the region properties of the objects
%                   rows,cols       : dimensions of the image
%
% output data:      NosplittedCells : those objects that were not splitted (small holes or one hole)
%                   splittedCells   : objects that were split
%                   NoHoleCells     : solid objects that did not contain holes
%                   These are presented as 3 separate 2D sets with the labelled objects
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


    
    
    if ~exist('rows') [rows,cols]           =  size(largeEndoCells); end
    if ~exist('statsObjects')
        %------- extracts the characteristics of the Remaining objects
        statsObjects                            =  regionprops(largeEndoCells, 'Area','EquivDiameter','Perimeter','Eccentricity','Centroid','BoundingBox','filledImage','filledarea','Image','EulerNumber');
    end
    %------- cycle over the large cells to determine their equivalente ellipse but first determine whether
    %------- A - they have several holes and thus need to be split into several cells or
    %------- B - they can be combined into a single new cell that was split before
    splittedCells(rows,cols)                =  0;
    NosplittedCells(rows,cols)              =  0;
    NoHoleCells(rows,cols)                  =  0;
    minAreaObject                           =  25;
    numLargeCells                           =  size(statsObjects,1);

    for counterObjs= 1:numLargeCells
        %------- if it has at least one hole split it
        %------- else, has no holes, then try to merge with neighbouring cells,
        %------- those that are within 2* its equivalent diameter
        minRow                              =  max(1,       ceil  (statsObjects(counterObjs).BoundingBox(2))           );
        maxRow                              =  min(rows,    floor (statsObjects(counterObjs).BoundingBox(4))+minRow -1 );
        minCol                              =  max(1   ,    ceil  (statsObjects(counterObjs).BoundingBox(1))           );
        maxCol                              =  min(cols,    floor (statsObjects(counterObjs).BoundingBox(3))+minCol -1 );
        %figure(3);surfdat(largeEndoCells==counterObjs);drawnow;
        %----- process only objects with holes
        if (statsObjects(counterObjs).EulerNumber<0)
            %----- areas are different if there is at least one hole in the object
            %----- label difference between image and the filled image, i.e. the HOLES
            [holesInData,numHolesInData]    = bwlabel(statsObjects(counterObjs).FilledImage-statsObjects(counterObjs).Image);
            if (numHolesInData>1)
                %------ more than one hole implies  linked objects that need to be splitted
                areaOfHoles                 = regionprops(holesInData,'Area');
                %------- remove holes that are smaller than a minimum Area **** 17 **** pixels
                indexToObjects              = find([areaOfHoles.Area] >=  minAreaObject);
                bigHolesOnly                = ismember(holesInData, indexToObjects);
                if length(indexToObjects)>1
                    %            [distToHoles,correspondingHoles]=  (bwdist(holesInData));
                    %[distToHoles,correspondingHoles]        =  (bwdist(bigHolesOnly));
                    [distToHoles,distBetHoles]              = regionGrowing(bigHolesOnly,1-statsObjects(counterObjs).FilledImage);
                    if any(distBetHoles(:,1)>3)

                        %------ avoid splitting between very close holes
                        distToHoles(distToHoles==0)             = max(distToHoles(:));
                        boundariesBetweenHoles                  =  watershed(distToHoles);
                        %------ close regions where distance between holes is less than 10
                        if any(distBetHoles(:,1)<10)
                            for k=1:size(distBetHoles)
                                if (distBetHoles(k,1)<10)&(sum(distBetHoles(k,4:5),2)<199)
                                    %------ last criterion, if holes are close, they need to have areas larger than 100 pixels each 
                                    boundariesBetweenHoles= boundariesBetweenHoles+uint8(imclose(ismember(boundariesBetweenHoles,[distBetHoles(k,2) distBetHoles(k,3) ]),[1 1 1;1 1 1; 1 1 1]));
                                end
                            end
                        end
                        [rBound,cBound]                         =  size(boundariesBetweenHoles);
                        %----- the boundaries may be too thin to separate the objects due to 4- 8- connectivity
                        %----- thicken the boundary
                        boundaryThick                           =  expandu(reduceu(boundariesBetweenHoles== 0)>0);
                        [SplittedObjects,numSplitedObjects]     =  bwlabel((1-boundaryThick(1:rBound,1:cBound)).*statsObjects(counterObjs).Image);

                        %                        [SplittedObjects,numSplitedObjects]=  bwlabel((boundariesBetweenHoles>0).*statsObjects(counterObjs).Image);
                        areaSplitedObjects=  regionprops(SplittedObjects,'Area');
                        indexToObjects  =   find([areaSplitedObjects.Area] >=  minAreaObject);
                        newObject=  ismember(SplittedObjects,indexToObjects);
                        splittedCells(minRow:maxRow,minCol:maxCol)=  newObject+splittedCells(minRow:maxRow,minCol:maxCol);
                        %splittedCells(minRow:maxRow,minCol:maxCol)=  SplittedObjects;
                        %                sC=  sC+1;
                    else
                        %------- holes are close to each other, leave as it is
                        NosplittedCells(minRow:maxRow,minCol:maxCol)=  statsObjects(counterObjs).Image+NosplittedCells(minRow:maxRow,minCol:maxCol);
                    end
                else
                    %------- several small holes, leave as it is
                    NosplittedCells(minRow:maxRow,minCol:maxCol)=  statsObjects(counterObjs).Image+NosplittedCells(minRow:maxRow,minCol:maxCol);
                    %               nSC=  nSC+1;
                end
            else
                %------- just one hole, leave as it is
                NosplittedCells(minRow:maxRow,minCol:maxCol)=  statsObjects(counterObjs).Image+NosplittedCells(minRow:maxRow,minCol:maxCol);
                %          nSC=  nSC+1;
            end
        else
            %-------- no holes, try to link with other neighbours in case is was a splitted cell
            NoHoleCells(minRow:maxRow,minCol:maxCol)=  NoHoleCells(minRow:maxRow,minCol:maxCol)+statsObjects(counterObjs).Image;
            %     nHC=  nHC+1;
        end
    end
    %        t1=  toc;tic;
    %        [sC nSC nHC numLargeCells]
    %        [cellsToLink,numCellsToLink]=  bwlabel(NoHoleCells);
    %        statsObjectsa=   regionprops(cellsToLink, 'Area','MajorAxisLength','Centroid');
    %        numLinked=  0;
    %index500   =  find(([statsObjectsa.Area]>500));
    %------ beginning of linking loop --------------------------------------------------------
    %         for counterObjs=  1:numCellsToLink-1
    %             %for counterObjs=  1:size(index500,2)-1
    %             %------ only link objects that are bigger than ... 100 pixels
    %             if statsObjectsa(counterObjs).Area>700
    %                 obj1=  (cellsToLink== counterObjs);
    %                 for counterObjs2=  counterObjs+1:numCellsToLink
    %                     %for counterObjs2=  1:size(index500,2)-1
    %                     %   index500(2:end)
    %                     if statsObjectsa(counterObjs2).Area>100
    %                         obj2=  (cellsToLink==  counterObjs2);
    %                         distBetCells=  sqrt(sum(diff([statsObjectsa(counterObjs).Centroid;statsObjectsa(counterObjs2).Centroid]).^2));
    %                         %------- only try to link if the cells are relatively close to each other
    %                         if (distBetCells<(0.5*statsObjectsa(counterObjs).MajorAxisLength+0.5*statsObjectsa(counterObjs2).MajorAxisLength))
    %                             [a,par,errorFit1,errorArc1] =  fitellipse(obj1);
    %                             [a,par,errorFit2,errorArc2] =  fitellipse(obj2);
    %                             [a,par,errorFit3,errorArc3] =  fitellipse(obj1+obj2);
    %                             if (errorArc3>errorArc2)&(errorArc3>errorArc1)
    %                                 if errorFit3<(0.5*errorFit1+0.5*errorFit2)
    %                                     numLinked= numLinked+1;
    %                                     params= a;
    %                                     %------ link both objects into one
    %                                     t =  linspace(0,pi*2,400);
    %                                     y =  params(4) * sin(t);
    %                                     x =  params(3) * cos(t);
    %                                     nx =  x*cos(params(5))-y*sin(params(5)) + params(1);
    %                                     ny =  x*sin(params(5))+y*cos(params(5)) + params(2);
    %                                     %------ the following loop  will draw a line over the two objects linking them
    %                                     for kkk= 1:400
    %                                         pos1= round(nx(kkk));
    %                                         pos2= round(ny(kkk));
    %                                         if (pos1>0)
    %                                             if (pos2>0)
    %                                                 cellsToLink(pos1,pos2)= counterObjs;
    %                                             end
    %                                         end
    %                                     end
    %                                     %------ this will change the values of the second objects to match the first
    %                                     cellsToLink(cellsToLink== counterObjs2)= counterObjs;
    %                                 end
    %                             end
    %                         end
    %                     end
    %                 end
    %             end
    %         end
    %
    %         t2= toc;tic;
    %        [numLinked numCellsToLink]

    %------ end of linking loop ----------------------------------------------------------------------
    %[finalCells,numNewObjects]= bwlabel((NosplittedCells+splittedCells+cellsToLink)>0);