% _______________________________________________________________________
%
% procosine.m (march, 8th 2016)
% _______________________________________________________________________


function PBRF = procosine(x,wl,Thetas,data)

% ***********************************************************************
% Jay, S., Bendoula, R., Hadoux, X., Féret, J.B. & Gorretta, N. (2016), A
% physically-based model for retrieving foliar biochemistry and leaf
% orientation using close-range imaging spectroscopy, Remote Sensing of 
% Environment, 177:220-236.
%
% Féret, J.B., François, C., Asner, G.P., Gitelson, A.A., Martin, R.E., 
% Bidel, L.P.R., Ustin, S.L., Le Maire, G. & Jacquemoud, S. (2008), PROSPECT-4 
% and 5: Advances in the Leaf Optical Properties Model Separating 
% Photosynthetic Pigments, Remote Sensing of Environment, 112:3030-3043.
% ***********************************************************************

N=x(1);
Cab=x(2);
Ccx=x(3);
Cbp=x(4);
Cw=x(5);
Cm=x(6);
Thetai=x(7);
Bspec=x(8);

RT = prospect_5B(N,Cab,Ccx,Cbp,Cw,Cm,data) ;

PBRF = cosine([Thetai Bspec],wl,RT(:,2),Thetas);

