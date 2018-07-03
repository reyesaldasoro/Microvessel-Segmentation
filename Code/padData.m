function dataPadded=padData (qtdata,numPadPixels,dimsToPad,padWith)
%function dataPadded=padData (qtdata,numPadPixels,dimsToPad,padWith)
%-----------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro-------------
%------             PHD     the University of Warwick-------------
%------  Supervisor :   Abhir Bhalerao    ------------------------
%------  5 March 2002    -----------------------------------------
%-----------------------------------------------------------------
%------ input  :  The data to be padded, default is to pad   -----
%------           with SAME values on the edges              -----
%------           but it can be changed with padWith to zero -----
%------           numPadPixels determine padding area        -----
%------ output :  the padded data
%-----------------------------------------------------------------
%----------------------------------------------------
%------ For a description and explanation please refer to:
%------ http://www.dcs.warwick.ac.uk/~creyes/m-vts --
%----------------------------------------------------
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

if ~exist('padWith') padWith=1; end
%-----dimensions will determine if it is a 1D, 2D, 3D, 4D ... padding
if (~exist('dimsToPad'))|(isempty(dimsToPad))
    [rows,cols,levs,numFeats]=size(qtdata);
else
    dimsToPad=[dimsToPad ones(1,4)]; %pad with ones in case it is too short
    rows=dimsToPad(1);
    cols=dimsToPad(2);
    levs=dimsToPad(3);
    numFeats=dimsToPad(4);
end

if rows>1    %---- first pad the rows
    qtdata=[padWith*repmat(qtdata(1,:,:,:),[numPadPixels 1 1]); qtdata;  padWith*repmat(qtdata(end,:,:,:),[numPadPixels 1 1])];
end
if cols>1    %---- then pad the cols
    qtdata=[padWith*repmat(qtdata(:,1,:,:),[1 numPadPixels 1]) qtdata  padWith*repmat(qtdata(:,end,:,:),[1 numPadPixels 1])];
end
dataPadded=qtdata;
if levs>1
    [rows,cols,levs,numFeats]=size(qtdata);
        qtdata3(:,:,numPadPixels+1:numPadPixels+levs,:)=qtdata;
        %qtdata3(:,:,numPadPixels+1:levs-numPadPixels,:)=qtdata;
        qtdata3(:,:,1:numPadPixels,:)=padWith*repmat(qtdata(:,:,1,:),[1 1 numPadPixels]);
        qtdata3(:,:,numPadPixels+levs+1:2*numPadPixels+levs,:)=padWith*repmat(qtdata(:,:,end,:),[1 1 numPadPixels]);
        %qtdata3(:,:,levs+1-numPadPixels:levs,:)=repmat(qtdata(:,:,end,:),[1 1 numPadPixels]);
        dataPadded=qtdata3;
end   
% if numFeats>1
%         for cFeats=1:numFeats
%             qtdata4(:,:,numPadPixels+1:numPadPixels+levs,cFeats)=qtdata3(:,:,:,cFeats); 
%             qtdata4(:,:,1:numPadPixels,cFeats)=repmat(qtdata3(:,:,1,cFeats),[1 1 numPadPixels]); 
%             qtdata4(:,:,numPadPixels+levs+1:2*numPadPixels+levs,cFeats)=repmat(qtdata3(:,:,end,cFeats),[1 1 numPadPixels]);
%         end
%         dataPadded=qtdata4;
% end

%         %---- convert into rows
%         for cFeats=1:numFeats qtdataCol(:,:,cFeats)=vol2col(qtdata4(:,:,:,cFeats),xSizeButterf*[1 1 1]); end
%         dataClassCol=vol2col(dataClass3,xSizeButterf*[1 1 1],boundary);
%         orientationTheta=orientationTheta(:);
%         orientationTheta(boundaryCol==0)=[];
%     else
%         %---- convert into rows
%         for cFeats=1:numFeats qtdataCol(:,:,cFeats)=im2colRed(qtdata3(:,:,cFeats),xSizeButterf*[1 1],boundary); end; 
%         dataClassCol=im2colRed(dataClass3,xSizeButterf*[1 1],boundary);
%     end
%     
%     
%         dataClass3=[repmat(dataClass(1,:,:,:),[numPadPixels 1 1]); dataClass;  repmat(dataClass(end,:,:,:),[numPadPixels 1 1])];
%     dataClass3=[repmat(dataClass3(:,1,:,:),[1 numPadPixels 1]) dataClass3  repmat(dataClass3(:,end,:,:),[1 numPadPixels 1])];
