function [h_hue_sat,h_hue_val,dataHue2,dataSaturation2,dataValue2]=colourHist2(dataHSV,sizeHue,sizeSaturation,sizeValue)
%function [h_hue_sat,h_hue_val,h_hue_sat_val,profMaxSat]=colourHist2(dataHSV,sizeHue,sizeSaturation,sizeValue)
%
%-------- this function obtains the joint Hue-Saturation and Hue-Value histograms together
%-------- with the Hue-Saturation-Value image to analyse the colour in images
%-------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------
%------             Postdoc  Sheffield University                           ----------
%------             http://tumour-microcirculation.group.shef.ac.uk         ----------
%------  29 May 2008   ---------------------------
%----------------------------------------------------
% input data:       dataHSV: an image with cells stained by immunohistochemsitry or any other colour image
%                   [sizeHue, sizeSaturation] : sizes of the histogram
% output data:      Corresponding Joint Histograms
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


%------ no input data is received, error -------------------------
%------ at least 2 parameters are required
if nargin <1;     help colourHist; h_hue_sat=[]; return;  end;
if ~exist('sizeHue','var');
    sizeHue         = 32;
    sizeSaturation  = 32;
    sizeValue       = 32;
end
if ~exist('sizeSaturation','var'); sizeSaturation  = ceil(sizeHue/2);                          end
if ~exist('sizeValue','var');      sizeValue       = ceil(sizeHue/2);                          end

%----- set H,S,V in different matrices
%----- if data has been previously quantised, uncomment following lines
dataHue                                     = dataHSV(:,:,1);
dataSaturation                              = dataHSV(:,:,2);
dataValue                                   = dataHSV(:,:,3);


[y,x]=hist(dataHue(:),(0:255)/255);
%if (numel(find(y)))~=sizeHue

    %%----- If data is not quantised, uncomment following lines
    dataHue             = 1+(sizeHue-1)*        quanti_r(dataHSV(:,:,1),log2(sizeHue));
    dataSaturation      = 1+(sizeSaturation-1)* quanti_r(dataHSV(:,:,2),log2(sizeSaturation));
    if max(max(dataHSV(:,:,3)))>1
        dataHSV(:,:,3)   = dataHSV(:,:,3)/255;
    end
    dataValue           = 1+(sizeValue-1)*      quanti_r(dataHSV(:,:,3),log2(sizeValue));


%else
%    dataHue             =  dataHue(:);
%    dataSaturation      =  dataSaturation(:);
%    dataValue           =  dataValue(:);

%end
dataHue2                                    = dataHue;
dataSaturation2                             = dataSaturation;
dataValue2                                  = dataValue;

chromaticity3D(sizeSaturation,sizeHue,sizeValue)=0;
%----- loop over the levels determined by the linspace xx
kk3=(1:sizeValue) ;
for k=1:sizeHue
    tempHue=find(dataHue==k);
    if ~isempty(tempHue)
        for k2=1:sizeSaturation
            %            tempSat=find(dataSaturation(tempHue)==k2);
            tempVal=(dataValue(tempHue));
            tempSat=tempVal(dataSaturation(tempHue)==k2);
            if ~isempty(tempSat)
                chromaticity3D(k2,k,kk3)=histc(tempSat,kk3);
            end
        end
        dataHue(tempHue)=[];
        dataSaturation(tempHue)=[];
        dataValue(tempHue)=[];
    end

end

%%
% dataHue3             =  dataHue(:);
% dataSaturation3      =  dataSaturation(:);
% dataValue3           =  dataValue(:);
% %%
% for k1=1:sizeHue
%     tempHue{k1} =find(dataHue3==k1);
% end
% for k1=1:sizeSaturation
%     tempSat{k1} =find(dataHue3==k1);
% end
% for k1=1:sizeValue
%     tempVal{k1} =find(dataHue3==k1);
% end
% %%
% 
% for k1=1:sizeHue
%     for k2=1:sizeSaturation
%         for k3=1:sizeValue
%             tempChrom =intersect(intersect(tempHue{k1},tempSat{k2}),tempVal{k3});
%             chromaticity3D(k2,k1,k3)=numel(tempChrom);
%         end
%     end
% end
% 
% for k1=1:sizeHue
%     for k2=1:sizeSaturation
%         for k3=1:sizeValue
%             tempChrom =numel(find((dataHue3==k1)&(dataSaturation3==k2)&(dataValue3==k3)));
%             chromaticity3D(k2,k1,k3)=tempChrom;
%         end
%     end
% end








h_hue_sat=sum(chromaticity3D,3);
h_hue_val=chromaticity3D;
