% _______________________________________________________________________
%
% cosine.m (march, 8th 2016)
% _______________________________________________________________________


function PBRF = cosine(x,wl,DHR,Thetas)

% ***********************************************************************
% Jay, S., Bendoula, R., Hadoux, X., Féret, J.B. & Gorretta, N. (2016), A
% physically-based model for retrieving foliar biochemistry and leaf
% orientation using close-range imaging spectroscopy, Remote Sensing of 
% Environment, 177:220-236.
% ***********************************************************************

Thetai=x(1);
Bspec=x(2);

PBRF = cosd(Thetai)/cosd(Thetas)*(DHR+Bspec);
