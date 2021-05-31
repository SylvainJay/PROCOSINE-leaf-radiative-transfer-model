% _______________________________________________________________________
%
% pro3cosine.m (july, 21th 2020)
% _______________________________________________________________________


function PBRF = pro3cosine(x,wl,Thetas,Coef_Spe)

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

% N=x(1);
% Cab=x(2);
% Cbp=x(3);
% Cw=x(4);
% Cm=x(5);
% Thetai=x(6);
% Bspec=x(7);

N=1.45;
Cab=x(1);
Cbp=0;
Cw=0.008;
Cm=0.005;
Thetai=x(2);
Bspec=x(3);

Ktot = (Coef_Spe(:,4)*Cab + Coef_Spe(:,5)*Cw ...
        + Coef_Spe(:,3)*Cm + Coef_Spe(:,6)*Cbp)./N;
RT = noyau(Coef_Spe(:,2),N,Ktot,tav(59.*pi/180,Coef_Spe(:,2))); % reflectance et transmittance des feuilles

PBRF = cosine([Thetai Bspec],wl,RT(:,2),Thetas);

