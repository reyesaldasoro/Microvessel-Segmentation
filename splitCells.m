function [splittedCellsLab,numNewObjects]=splitCells(dataIn,minAreaObject)
%function [splittedCellsLab,numNewObjects]=splitCells(dataIn,minAreaObject)
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


%------- function to split a number of cells if they have more than one inner region
[rows,cols,levs]=size(dataIn);
splittedCells=zeros(rows,cols);
maxData=max(dataIn(:));
minData=min(dataIn(:));

if ~exist('minAreaObject')        minAreaObject=17; end
if (isa(dataIn,'logical'))|((minData==0)&(maxData==1))
    %----- this is the case of a binary image, labelling is required
    [dataLabelled,numObjects]=bwlabel(dataIn,8);
else
    %----- in this case the image has more than 2 levels, do not label
    numObjects=maxData;
end

for counterObjs=1:numObjects
    %--------- ??????? filter to remove small dots??????
    %filtG=gaussF(3,3,1);
    %dataLabelled=conv2(double(handles.labelledData==1363),filtG,'same')>0.1;
    %--------- fill the current object to see if it has more than one "hole"
    currentObject=(  dataLabelled==counterObjs );
    rowProj=find(max(currentObject'));
    colProj=find(max(currentObject));
    minRow=max(rowProj(1)-1,1);
    maxRow=min(rowProj(end)+1,rows);
    minCol=max(colProj(1)-1,1);
    maxCol=min(colProj(end)+1,cols);
    currentObject2=currentObject(minRow:maxRow,minCol:maxCol);
    dataFilled=(imfill(currentObject2,'holes'));
    [holesInData,numHolesInData]=bwlabel(dataFilled-currentObject2);
    if numHolesInData>1
        areaOfHoles=regionprops(holesInData,'Area');
        %------- remove holes that are smaller than **** 20 **** pixels

        indexToObjects  = find([areaOfHoles.Area] >=minAreaObject);
        bigHolesOnly     = ismember(holesInData, indexToObjects);

        if length(indexToObjects)>1
            %------ more than one hole implies  linked objects that need to be splitted

            [distToHoles,correspondingHoles]=(bwdist(bigHolesOnly));
            boundariesBetweenHoles=watershed(distToHoles);
            [SplittedObjects,numSplitedObjects]=bwlabel((boundariesBetweenHoles>0).*currentObject2);
            areaSplitedObjects=regionprops(SplittedObjects,'Area');
            indexToObjects  = find([areaSplitedObjects.Area] >=minAreaObject);
            newObject=ismember(SplittedObjects,indexToObjects);
            splittedCells(minRow:maxRow,minCol:maxCol)=newObject;
         %           figure(11);subplot(211);surfdat(currentObject2);subplot(212);surfdat(bwlabel(newObject))
%qqqq=0;
        else
            splittedCells(minRow:maxRow,minCol:maxCol)=currentObject2;

        end

    else
        splittedCells(minRow:maxRow,minCol:maxCol)=currentObject2;
    end

end

[splittedCellsLab,numNewObjects]=bwlabel(splittedCells);