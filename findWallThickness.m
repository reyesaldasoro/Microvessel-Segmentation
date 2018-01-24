function  [avThickness,stdThickness,propLumen,minThickness,maxThickness]=findWallThickness(dataIn)


%-------- regular size check and determination of first zoom region
[rows,cols]         = size(dataIn);
%surfdat(dataIn)
%-------- first crop image to the region where the object is, this saves
%-------- computing other properties for a big image
if isa(dataIn,'logical')
    dataIn              = double(dataIn>0);
    %    r1                  = max(1,   floor (tempProps.BoundingBox(2)));
    %    r2                  = min(rows,ceil  (tempProps.BoundingBox(2)+tempProps.BoundingBox(4)));
    %    c1                  = max(1,   floor (tempProps.BoundingBox(1)));
    %    c2                  = min(cols,ceil  (tempProps.BoundingBox(1)+tempProps.BoundingBox(3)));
    %    dataInR              = (dataIn(r1:r2,c1:c2));
    [avThickness,stdThickness,propLumen,minThickness,maxThickness]=findWallThickness(dataIn);
else
    
    if max(dataIn(:))>1
        tempProps1          = regionprops(dataIn,'boundingbox');
        numObjects          = length(tempProps1);
        avThickness=zeros(numObjects,1);minThickness=zeros(numObjects,1);
        stdThickness=zeros(numObjects,1);maxThickness=zeros(numObjects,1);
        propLumen=zeros(numObjects,1);
        for counterObj=1:numObjects
            %counterObj
           
            r1                  = max(1,   floor (tempProps1(counterObj).BoundingBox(2)));
            r2                  = min(rows,ceil  (tempProps1(counterObj).BoundingBox(2)+tempProps1(counterObj).BoundingBox(4)));
            c1                  = max(1,   floor (tempProps1(counterObj).BoundingBox(1)));
            c2                  = min(cols,ceil  (tempProps1(counterObj).BoundingBox(1)+tempProps1(counterObj).BoundingBox(3)));
            
            [avThickness(counterObj),stdThickness(counterObj),propLumen(counterObj),minThickness(counterObj),maxThickness(counterObj)]=findWallThickness(double(dataIn(r1:r2,c1:c2)==counterObj));
        end
    else

        %------- calculate now other properties that will be required
        tempProps           = regionprops(dataIn,'FilledArea','Image','FilledImage');
        
        %------- morphologically THIN the image and remove a few spurious lines
        %------- only for relatively big objects
        centreLine          = bwmorph(tempProps.Image,'thin','inf');
        
        %------- find inner holes or lumen of the object and measure the area
        [holesInData]    = bwlabeln(tempProps.FilledImage-tempProps.Image);

        %------- generate distance map from the outside to grow into the object, in
        %------- that way, the map represents the distance from the outer boundary
        distMap             = (bwdist(1-padData(tempProps.FilledImage,1,[],0),'euclidean'));
        %------- remove the holes from the distance map, not sure if necessary
        distMap             = distMap(2:end-1,2:end-1).*tempProps.Image;

 
        if sum(centreLine(:))>20
            centreLine      = bwmorph(centreLine,'spur',5);
        end
        %------ the line is times 2 as centre line has two sides
        weightedCentreLine  = (2*(centreLine.*distMap));
        %------ get positions of the valid pixels and turn into a vector
        positionsCentreLine = (find(weightedCentreLine));
        vectorThickness     = weightedCentreLine(positionsCentreLine);
        %------ finally calculate all parameters
        avThickness         = mean (vectorThickness);
        stdThickness        = std  (vectorThickness);
        minThickness        = min  (vectorThickness);
        maxThickness        = max  (vectorThickness);
        propLumen           = sum(holesInData(:)>0)/tempProps.FilledArea;
       % subplot(131);surfdat(dataIn);subplot(132);surfdat((1-centreLine).*distMap);subplot(133);surfdat(weightedCentreLine);

    end
end


