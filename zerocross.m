function bor = zerocross(imag)
%function bor = zerocross(imag)

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



% get zero crossings of a certain input can be 2D or 1D

% if  min(imag(:))>= 0;
%    msgbox('No negative values in data','Zerocross warning','warn');
%    bor = [];
%    return;
% end;

[lins,cols,levels] = size(imag);
delta = 0.00001;

if ~isa(imag,'logical')
    imag=sign(imag);
end
%------ revise the case ------
%------    1 1D use plot for either line or column data
%------    2 1D but not line or column, stored in various z
%------    3 2D use surfdat_r.m function
%------    4 3D use just the base of cube

if ((cols==1|lins==1)&levels==1)  %- 1D over main plane
    if cols>= lins
        yy = [0 imag(1:cols-1)];
    else
        yy = [0 imag(1:lins-1)']';
    end;
    bor = abs((sign(imag+delta)-sign(yy+delta)));
elseif (cols==1&lins==1&levels~= 1)%- 1D over other plane
    imag = permute(imag,[2 3 1]);
    yy = [0 imag(1:cols-1)];
    bor = (sign(imag+delta)-sign(yy+delta));
elseif (lins~= 1&cols~= 1&levels==1) %- 2D over main plane
    %------ only 1 degree neighbourhood considered---------
    %    imag = imag+delta;
    %    yy5 = [zeros(1,cols);  imag(1:lins-1,1:cols)];              %|d
    %    yy6 = [imag(2:lins,1:cols);zeros(1,cols)];                  %u|
    %    yy7 = [ zeros(lins,1) imag(1:lins,1:cols-1)];               %-r
    %    yy8 = [ imag(1:lins,2:cols) zeros(lins,1)];                 %l-
    %    bor5 = fix(delta+(sign(imag)-sign(yy5))/2);
    %    bor6 = fix(delta+(sign(imag)-sign(yy6))/2);
    %    bor7 = fix(delta+(sign(imag)-sign(yy7))/2);
    %    bor8 = fix(delta+(sign(imag)-sign(yy8))/2);
    %    bor = sign(bor5+bor6+bor7+bor8);
    diffVer         = diff(imag,1,1);zerCols = zeros(1,cols);
    diffHor         = diff(imag,1,2);zerRows = zeros(lins,1);
    qq1             = [zerCols;(diffVer)>0];
    qq2             = [(diffVer)<0;zerCols];
    qq3             = [ (diffHor)<0 zerRows ];
    qq4             = [ zerRows (diffHor)>0 ];
    bor             = qq1|qq2|qq3|qq4;
elseif(lins~= 1&cols~= 1&levels~= 1) %- 3D
    yy5             = [zeros(1,cols,levels);  imag(1:lins-1,1:cols,:)];     %|d
    yy6             = [imag(2:lins,1:cols,:);zeros(1,cols,levels)];                  %u|
    yy7             = [ zeros(lins,1,levels) imag(1:lins,1:cols-1,:)];               %-r
    yy8             = [ imag(1:lins,2:cols,:) zeros(lins,1,levels)];                 %l-
    bor5            = fix(delta+(sign(imag)-sign(yy5))/2);
    bor6            = fix(delta+(sign(imag)-sign(yy6))/2);
    bor7            = fix(delta+(sign(imag)-sign(yy7))/2);
    bor8            = fix(delta+(sign(imag)-sign(yy8))/2);
    bor             = sign(bor5+bor6+bor7+bor8);
end;

%------ 2nd degree neighbourhood ----------------------------
%   yy1 = [zeros(1,cols); zeros(lins-1,1) imag(1:lins-1,1:cols-1)]; %\d
%   yy2 = [zeros(1,cols); imag(1:lins-1,2:cols) zeros(lins-1,1)];   %/d
%   yy3 = [imag(2:lins,2:cols) zeros(lins-1,1); zeros(1,cols)];     %u\
%   yy4 = [zeros(lins-1,1) imag(2:lins,1:cols-1); zeros(1,cols)];   %u/
%   bor1 = abs(sign(imag+delta)-sign(yy1+delta));
%   bor2 = abs(sign(imag+delta)-sign(yy2+delta));
%   bor3 = abs(sign(imag+delta)-sign(yy3+delta));
%   bor4 = abs(sign(imag+delta)-sign(yy4+delta));
%   bor = sign(bor1+bor2+bor3+bor4+bor5+bor6+bor7+bor8);
