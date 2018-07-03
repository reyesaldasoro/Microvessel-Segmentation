function [dataOut,distBetweenHoles] = regionGrowing (initialSeeds,backgroundMask)

    %-------- regular size check and determination of first zoom region
    [rows,cols]=size(initialSeeds);

    if ~exist('backgroundMask','var')
        backgroundMask=zeros(rows,cols);
    end

    if isa(initialSeeds,'logical')
        initialSeeds                = bwlabel(initialSeeds);
    end

    numObjects=max(initialSeeds(:));

    if numObjects>1
        
        if numObjects>100
            try
            dataOut2=zeros(rows,cols,100);
            catch
                qqq=1;
            end
            numIterations                   = ceil(numObjects/100);
            dataOut3=zeros(rows,cols,numIterations);
            for countIteration = 1: numIterations
                disp(strcat('numiterations:',num2str(countIteration)));
                counterLowLevel             = 100*(countIteration-1);
                for counterObjs=counterLowLevel+1:min(100,numObjects-counterLowLevel)
                    [dataOut2(:,:,counterObjs)]  = regionGrowing (initialSeeds==counterObjs,backgroundMask);
                end
                if nargout>1
                    distBetweenHoles            = [];
                    counterLowLevel             = 100*(countIteration-1);
                    for counterHoles=counterLowLevel+1:min(100,numObjects-counterLowLevel)-1
                        for counterHoles2=counterHoles+1:min(100,numObjects-counterLowLevel)
                            tempDist            = dataOut2(:,:,counterHoles).*(initialSeeds==(counterHoles2));
                            try
                                %[distance Hole1  Hole2  SizeHole1 SizeHole2]
                                distBetweenHoles= [distBetweenHoles;[min(tempDist(find(tempDist))) counterHoles  counterHoles2 sum(find(initialSeeds==counterHoles)>0) sum(find(initialSeeds==counterHoles2)>0)]];
                            catch
                                ttt=1;
                            end

                        end
                    end
                end
                dataOut3(:,:,countIteration)    = min(dataOut2,[],3);


            end

            dataOut=min(dataOut3,[],3);
            
            
            
        else
            dataOut=zeros(rows,cols,numObjects);
            for counterObjs=1:numObjects
                [dataOut(:,:,counterObjs)]  = regionGrowing (initialSeeds==counterObjs,backgroundMask);

            end
            if nargout>1
                distBetweenHoles            = [];
                for counterHoles=1:(numObjects-1)
                    for counterHoles2=counterHoles+1:numObjects
                        tempDist            = dataOut(:,:,counterHoles).*(initialSeeds==(counterHoles2));
                        try
                            %[distance Hole1  Hole2  SizeHole1 SizeHole2]
                            distBetweenHoles= [distBetweenHoles;[min(tempDist(find(tempDist))) counterHoles  counterHoles2 sum(find(initialSeeds==counterHoles)>0) sum(find(initialSeeds==counterHoles2)>0)]];
                        catch
                            ttt=1;
                        end

                    end
                end
            end


            dataOut=min(dataOut,[],3);
        end
    else
        fourConnectedKernel         = strel    ([0 1 0;1 1 1; 0 1 0]);
        outerBoundary               = -initialSeeds+imdilate (initialSeeds,fourConnectedKernel);
        seedsProps                  = find(initialSeeds);
        outerB_Props                = find(outerBoundary);
        %seedsProps                  = regionprops(initialSeeds, 'area','pixellist','pixelidxlist');
        %outerB_Props                = regionprops(outerBoundary,'area','pixellist','pixelidxlist');
        %------- define outer Boundary
        outerBLocation           = outerB_Props;
        seedLocation=find(initialSeeds(:));
        %%
        k=0;
        while( ~isempty(outerBLocation) )
            k=k+1;
            %----- some pixels are outside the image
            outerBLocation(outerBLocation<1)=[];
            outerBLocation(outerBLocation>(rows*cols))=[];
            % First, discard all pixels that belong to the background Mask
            try
                outerBLocation(backgroundMask(outerBLocation)>0)=[];
            catch

            end
            if ~isempty (outerBLocation)
                %------- growing process until it stops
                % update regions include new pixel on initialSeeds,  join the pixels that is closest to the current Seed
                initialSeeds(outerBLocation) = k;
%                figure(2);surfdat(initialSeeds);axis off;drawnow;
                %jj=length(FRA_1)+1
                %FRA_1(jj)=getframe;
                seedsProps=[seedsProps;outerBLocation];

                %----- add pixels added "outerBlocation" into the seedLocation
                seedLocation                                        = [seedLocation;outerBLocation];
                %----- find new neighbours, up, down, left and right
                newNeighboursUp                                     = outerBLocation-1;
                newNeighboursDown                                   = outerBLocation+1;
                newNeighboursRight                                  = outerBLocation+rows;
                newNeighboursLeft                                   = outerBLocation-rows;
                %----- find new neighbours, NE, NW, SE and SW
                newNeighboursNE                                     = newNeighboursRight-1;
                newNeighboursSE                                     = newNeighboursRight+1;
                newNeighboursNW                                     = newNeighboursLeft-1;
                newNeighboursSW                                     = newNeighboursLeft+1;
                
                
                
                %----- check that none of the neighbours are outside the image
%                 boundCheck                                          = mod(outerBLocation,rows);
%                  newNeighboursUp    (boundCheck==1)                  =[];
%                  newNeighboursDown  (boundCheck==0)                  =[];
%                  newNeighboursRight (outerBLocation>(rows*(cols-1))) =[];
%                  newNeighboursLeft  (outerBLocation<rows)            =[];
%                 newNeighboursNE (boundCheck==1)                  =[];
%                 newNeighboursSE (boundCheck==0)                  =[];
%                 newNeighboursNW (outerBLocation>(rows*(cols-1))) =[];
%                 newNeighboursSW (outerBLocation<rows)            =[];
              
                 
                 
                newNeighbours4Conn = [newNeighboursRight;newNeighboursLeft;newNeighboursUp;newNeighboursDown];
                newNeighbours8Conn = [newNeighboursNE;newNeighboursSE;newNeighboursNW;newNeighboursSW];
                newNeighbours = [newNeighbours4Conn;newNeighbours8Conn];

                newNeighbours(newNeighbours<1) = [];
                newNeighbours(newNeighbours>(rows*cols)) = [];
                newNeighbours(mod(newNeighbours,rows)==0) = [];
                newNeighbours(mod(newNeighbours,rows)==1) = [];
                %----- add new neighbours into the outerBLocation
                %outerBLocation= unique( [outerBLocation;newNeighboursRight;newNeighboursLeft;newNeighboursUp;newNeighboursDown]);
                outerBLocation= unique( [outerBLocation;newNeighbours]);
                %----- remove from OuterBLocation all the pixels that belong to seedLocation
                outerBLocation(ismember(outerBLocation,seedLocation))=[];
            end
        end% end of the loop
        dataOut=initialSeeds;

        %         if nargout>1
        %
        %             tempDist=dataOut(:,:,counterHoles).*(initialSeeds==(counterHoles+1));
        %             distBetweenHoles=[min(tempDist(find(tempDist)))];
        %
        %
        %         end

    end