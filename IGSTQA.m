% This code is written by: S. Alireza Golestaneh;
% Email: sgolest1@asu.edu;
% March 1 2017; for CVPRW 2018; Paper: Synthesized Texture Quality Assessment 
% via Multi-scale Spatial and Statistical Texture Attributes of Image and Gradient Magnitude Coefficients
% S. Alireza Golestaneh and Lina Karam

% @InProceedings{Golestaneh_2018_CVPR,
% author = {Alireza Golestaneh, S. and Karam, Lina J.},
% title = {Spatially-Varying Blur Detection Based on Multiscale Fused and Sorted Transform Coefficients of Gradient Magnitudes},
% booktitle = Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition Workshops},
% month = {July},
% year = {2018}
% }

% This code was tested on MATLAB2017a

%=====================================================================
%%
% Main Funcrion
function QualityScore = IGSTQA(Ref_Img,Synth_Img)
% Making Sure images are in grayscale domain
if length(size(Synth_Img))==3
    Synth_Img=rgb2gray(Synth_Img);
end
if length(size(Ref_Img))==3
    Ref_Img=double(rgb2gray(Ref_Img));
end

% Reference and Synthesized images
Ref_Img = double(Ref_Img);
Synth_Img = double(Synth_Img);

% Computing the gradient magnitudes (GM) for the Reference and Synthesized images
Synth_Img_GM = double(imgradient(Synth_Img));
Ref_Img_GM = double(imgradient(Ref_Img));

% Number of scales that are used for the wavelet
Scales = 4;

% Statistical Features in Vertical bands while using the spatial-domain image as an input --> Ref_F_V
% Statistical Features in Horizental bands while using the spatial-domain image as an input --> Ref_F_H
% Distance Features in Vertical bands while using the spatial-domain image as an input --> Ref_D_V
% Distance Features in Horizental bands while using the spatial-domain image as an input --> Ref_D_H
[Ref_F_V Ref_F_H  Ref_D_V Ref_D_H] = Computing_Features(Ref_Img,Scales);
[Dist_F_V  Dist_F_H  Dist_D_V Dist_D_H] = Computing_Features(Synth_Img,Scales);

% Statistical Features in Vertical bands while using the GM image as an input --> Ref_F_V_GM
% Statistical Features in Horizental bands while using the GM image as an input --> Ref_F_H_GM
% Distance Features in Vertical bands while using the GM image as an input --> Ref_D_V_GM
% Distance Features in Horizental bands while using the GM image as an input --> Ref_D_H_GM
[Ref_F_V_GM Ref_F_H_GM  Ref_D_V_GM Ref_D_H_GM] = Computing_Features(Ref_Img_GM,Scales);
[Dist_F_V_GM Dist_F_H_GM  Dist_D_V_GM Dist_D_H_GM] = Computing_Features(Synth_Img_GM,Scales);


%%
% Computing Statistical Features, Granularity, and Regularity for the spatial domain input
Statistical_F_Differences_Vertical =  mean(abs(Ref_F_V-Dist_F_V),1);
Statistical_F_Differences_Horizental =  mean(abs(Ref_F_H-Dist_F_H),1) ;
Statistical_F_Differences = (Statistical_F_Differences_Vertical + Statistical_F_Differences_Horizental)/2;
for Scale_Num = 1:Scales
    Granularity_Vertical(Scale_Num,1) = abs(mean(Ref_D_H{Scale_Num})-mean(Dist_D_H{Scale_Num}));
    Granularity_Horizental(Scale_Num,1) = abs(mean(Ref_D_V{Scale_Num})-mean(Dist_D_V{Scale_Num}));
    
    Regularity_Vertical(Scale_Num,1) = abs(std(Ref_D_H{Scale_Num})-std(Dist_D_H{Scale_Num}));
    Regularity_Horizental(Scale_Num,1) = abs(std(Ref_D_V{Scale_Num})-std(Dist_D_V{Scale_Num}));
end
Total_Granularity_Difference  = (Granularity_Vertical  + Granularity_Horizental)/2;
Total_Regularity_Difference   =  (Regularity_Vertical  + Regularity_Horizental)/2;
%%
% Computing Statistical Features, Granularity, and Regularity for the GM domain input
Statistical_F_Differences_Vertical_GM =  mean(abs(Ref_F_V_GM-Dist_F_V_GM),1);
Statistical_F_Differences_Horizental_GM =  mean(abs(Ref_F_H_GM-Dist_F_H_GM),1);
Statistical_F_Differences_GM = (Statistical_F_Differences_Vertical_GM +Statistical_F_Differences_Horizental_GM)/2;

for Scale_Num=1:Scales
    Granularity_Vertical_GM(Scale_Num,1) = abs(mean(Ref_D_H_GM{Scale_Num})-mean(Dist_D_H_GM{Scale_Num}));
    Granularity_Horizental_GM(Scale_Num,1) = abs(mean(Ref_D_V_GM{Scale_Num})-mean(Dist_D_V_GM{Scale_Num}));
    
    Regularity_Vertical_GM(Scale_Num,1) = abs(std(Ref_D_H_GM{Scale_Num})-std(Dist_D_H_GM{Scale_Num}));
    Regularity_Horizental_GM(Scale_Num,1) = abs(std(Ref_D_V_GM{Scale_Num})-std(Dist_D_V_GM{Scale_Num}));
end
Total_Granularity_Difference_GM  = (Granularity_Vertical_GM  + Granularity_Horizental_GM )/2;
Total_Regularity_Difference_GM =  (Regularity_Vertical_GM  + Regularity_Horizental_GM)/2;

%%
% Computing Quality Features for Spatial and GM domains
QualityFeatures_Spatial=log(1+(100.*[ Statistical_F_Differences   max(Total_Regularity_Difference',[],2,'omitnan')  max(Total_Granularity_Difference',[],2,'omitnan')   ]));
QualityFeatures_Spatial_GM =log(1+(100.*[ Statistical_F_Differences_GM   max(Total_Regularity_Difference_GM',[],2,'omitnan')  max(Total_Granularity_Difference_GM',[],2,'omitnan')   ]));

% Total Quality Score
QualityScore= mean(QualityFeatures_Spatial,2)+mean(QualityFeatures_Spatial_GM,2) ;
end
%%
% Other Required Functions:
function [STAT_L_O_2 STAT_H_O_2 PEAKS_L_O_2 PEAKS_H_O_2]=Computing_Features(A0_IMG,scales)

% [Stat_V_A Stat_H_A  Peak_V_A Peak_H_A] = [STAT_L_O_2 STAT_H_O_2 PEAKS_L_O_2 PEAKS_H_O_2]
% Statistical Features in Vertical bands while using the spatial-domain image as an input --> Ref_F_V
% Statistical Features in Horizental bands while using the spatial-domain image as an input --> Ref_F_H
% Distance Features in Vertical bands while using the spatial-domain image as an input --> Ref_D_V
% Distance Features in Horizental bands while using the spatial-domain image as an input --> Ref_D_H

% computing wavelets
[im_ll_O, im_lh_O, im_hl_O] = FWT_ATrou(A0_IMG , scales);
for k=1:scales
    
    % computing features and peaks for im_lh_O
    MATR_L_O_2 = abs(im_lh_O(:,:,k) );
    STAT_L_O_2(k,:) = [ std2(MATR_L_O_2) kurtosis(MATR_L_O_2(:)) skewness(MATR_L_O_2(:)) wentropy(MATR_L_O_2,'log energy') ];
    PEAKS_L_O_2{k,1} = findpeaks( (MATR_L_O_2(:)),'Threshold', max(mean2( (MATR_L_O_2))/1,0));
    
    % computing features and peaks for im_hl_O
    MATR_H_O_2  =       abs(im_hl_O(:,:,k) );
    STAT_H_O_2(k,:) = [  std2(MATR_H_O_2)   kurtosis(MATR_H_O_2(:)) skewness(MATR_H_O_2(:)) wentropy(MATR_H_O_2,'log energy')   ];
    PEAKS_H_O_2{k,1} = findpeaks( (MATR_H_O_2(:)),'Threshold', max( mean2( (MATR_H_O_2))/1,0));
end
end
%%
function [A,H,V] = FWT_ATrou(x,L);
% FWT_ATrou -- Fast Dyadic Wavelet Transform (periodized, orthogonal)
%  Usage:
%  dwt = FWT_ATrou(x,L)
%  Inputs:
%    x    	1-d signal; length(x) = 2^J = n
%    L    	Coarsest Level of V_0;  L << J
%  Outputs:
%    dwt   an n times J-L+1 matrix
%           giving the wavelet transform of x at all dyadic scales.
%  Description
%    To reconstruct use IWT_ATrou
%  See Also:
%    IWT_ATrou, MakeATrouFilter
%
%[lodyadf,dlodyadf,hidyadf,dhidyadf] = MakeATrouFilter('Spline',3);
lodyadf = [0.125 0.375 0.375 0.125];
hidyadf = [-2.0 2.0];
lamda = [1.50 1.12 1.03 1.01 1 1 1 1 1 1 1];
[n,m] = size(x);
D = L;

current_low = x;
shiftH = 0;
shiftV = 0;

for d=1:D
    % perform the row filtering first
    for i=1:1:n
        s=ShapeAsRow(current_low(i,:));
        % perform high frequency filtering for the row
        H(i,:,d) = iconv(hidyadf,s)./lamda(d);
        % do the shift, why?
        %for j = 1:2^(d)
        %  p = lshift(H(i,:,d));
        %  H(i,:,d) = p;
        %end
        % perform the low frequency filter for the row
        s2 = s;
        for j = 1:2^(d-1)
            s2 = lshift(s2);
        end
        A(i,:,d) = iconv(lodyadf,s2);
    end
    % perform the column filtering
    for i=1:1:m
        s=ShapeAsRow(current_low(:,i));
        % perform high frequency filtering for the column
        V(:,i,d) = iconv(hidyadf,s)'./lamda(d);
        % do the shift, why?
        % for j = 1:2^(d)
        %     p = lshift(V(:,i,d)');
        %     V(:,i,d) = p';
        % end
        % perform the low frequency filter for the row
        s2 = ShapeAsRow(A(:,i,d));
        for j = 1:2^(d-1)
            s2 = lshift(s2);
        end
        A(:,i,d) = iconv(lodyadf,s2)';
    end
    
    % shift A(:,:,d)
    for row=1:1:n
        for i=1:1:shiftH
            A(row,:,d)=lshift(A(row,:,d));
            H(row,:,d)=lshift(H(row,:,d));
        end
    end
    for column=1:1:m
        for i=1:1:shiftV
            A(:,column,d)=lshift(A(:,column,d)')';
            V(:,column,d)=lshift(V(:,column,d)')';
        end
    end
    shiftH=2^(d-1);
    shiftV=2^(d-1);
    current_low=A(:,:,d);
    f = zeros(1,2*length(lodyadf));
    f(1:2:2*length(lodyadf)-1) = lodyadf;
    f2 = zeros(1,2*length(hidyadf));
    f2(1:2:2*length(hidyadf)-1) = hidyadf;
    lodyadf = f;
    hidyadf = f2;
end
end
%%
function row = ShapeAsRow(sig)
% ShapeAsRow -- Make signal a row vector
%  Usage
%    row = ShapeAsRow(sig)
%  Inputs
%    sig     a row or column vector
%  Outputs
%    row     a row vector
%
%  See Also
%    ShapeLike
%
row = sig(:)';    

%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:43 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
end
%%
function y = iconv(f,x)
% iconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv(f,x)
%  Inputs
%    f   filter
%    x   1-d signal
%  Outputs
%    y   filtered result
%
%  Description
%    Filtering by periodic convolution of x with f
%
%  See Also
%    aconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%
	n = length(x);
	p = length(f);
	if p <= n,
	   xpadded = [x((n+1-p):n) x];
	else
	   z = zeros(1,p);
	   for i=1:p,
		   imod = 1 + rem(p*n -p + i-1,n);
		   z(i) = x(imod);
	   end
	   xpadded = [z x];
	end
	ypadded = filter(f,1,xpadded);
	y = ypadded((p+1):(n+p));
%
% Copyright (c) 1993. David L. Donoho
%     

%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
end
%%
function y = lshift(x)
% lshift -- Circular left shift of 1-d signal
%  Usage
%    l = lshift(x)
%  Inputs
%    x   1-d signal
%  Outputs
%    l   1-d signal 
%        l(i) = x(i+1) except l(n) = x(1)
%
	y = [ x( 2:length(x) ) x(1) ];
%
% Copyright (c) 1993. Iain M. Johnstone
%     

%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
end
%%