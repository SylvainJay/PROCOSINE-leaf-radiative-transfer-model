% _______________________________________________________________________
% main.m
% version 1 (March, 8th 2016)
% subroutines required: cosine.m, procosine.m
%                       prospect_5B.m, tav.m, dataSpec_P5B.m
% Matlab toolbox required: Optimization Toolbox  
% author: Sylvain Jay (sylvain.jay@fresnel.fr) 
% _______________________________________________________________________
% 
% This script allows inverting the PROSPECT+COSINE model based on a pseudo 
% bidirectional reflectance factor measurement and iterative optimization.
% In the following, PROSPECT-5 (Féret et al., 2008) is used to simulate 
% the leaf directional-hemispherical reflectance in the optical domain 
% (400-2500 nm) with 1 nm step. The matlab codes of PROSPECT-5 (and PROSPECT-4) 
% are available on the OPTICLEAF website at 
% http://teledetection.ipgp.jussieu.fr/prosail/.
% Inverting the combined PROCOSINE model allows the estimation of the
% following 8 parameters:
%
%       - N      = leaf structure parameter
%       - Cab    = chlorophyll a+b content in µg/cm²
%       - Car    = carotenoids content in µg/cm²
%       - Cbrown = brown pigments content in arbitrary units
%       - Cw     = equivalent water thickness in g/cm² or in cm
%       - Cm     = dry matter content in g/cm²
%       - Thetai = light incident angle in °
%       - Bspec  = specular parameter (unitless)
%
% Model inversion is performed using the trust-region reflective algorithm
% implemented in the matlab function lsqcurvefit.m%
% _______________________________________________________________________
% 
% Jay, S., Bendoula, R., Hadoux, X., Féret, J.B. & Gorretta, N. (2016), A
% physically-based model for retrieving foliar biochemistry and leaf
% orientation using close-range imaging spectroscopy, Remote Sensing of
% Environment, 177:220-236.
%
% Féret, J.B., François, C., Asner, G.P., Gitelson, A.A., Martin, R.E., 
% Bidel, L.P.R., Ustin, S.L., Le Maire, G., Jacquemoud, S. (2008), PROSPECT-4 
% and 5: Advances in the Leaf Optical Properties Model Separating 
% Photosynthetic Pigments, Remote Sensing of Environment, 112:3030-3043.
% _______________________________________________________________________


clear all
clc

%% Model simulation

dir *.txt
filename=input('file name (without extension): ','s'); % PROCOSINE parameters (leaf_parameters.txt)
eval(['load ',filename,'.txt']);
eval(['leaf_parameters=',filename,';']);
data= dataSpec_P5B;

Thetas = input('Illumination zenith angle Theta_s (in degree) ?');

rsim = procosine(leaf_parameters,data(:,1),Thetas,data);

figure,plot(data(:,1),rsim,'LineWidth',2.5,'Color',[0.75,0.75,0.99]);
xlabel('Wavelength (nm)','FontSize',11,'FontWeight','Bold')
ylabel('Pseudo Bidirectional Reflectance Factor','FontSize',12,'FontWeight','Bold')
title(['N = ',num2str(leaf_parameters(1),'%4.2f'),'  C_{ab} = ',num2str(leaf_parameters(2),'%4.1f'),...
    '  C_{cx} = ',num2str(leaf_parameters(3),'%4.2f'),'  C_{bp} = ',num2str(leaf_parameters(4),'%4.1f'),...
    '  C_{w} = ',num2str(leaf_parameters(5),'%7.4f'),'  C_{dm} = ',num2str(leaf_parameters(6),'%7.4f'),...
    '  \theta_{i} = ',num2str(leaf_parameters(7),'%4.1f'),'  b_{spec} = ',num2str(leaf_parameters(8),'%4.2f')],'FontWeight','Bold')
axis([400 2500 0 1])




%% Model inversion

dir *.mat
filename=input('file name (without extension): ','s'); % Measured pseudo-BRF data
eval(['load ',filename,'.mat']); % Load pseudo-BRF image, camera wavelength and Thetas
im_pbrf = double(im_pbrf)/10000;
rgb = [55 41 12];
figure,imshow(im_pbrf(:,:,rgb)*1.5);
rmes=squeeze(im_pbrf(25,64,:)); % Spectrum coordinates in the image to replace with the desired ones

data= dataSpec_P5B;
data = spline(data(:,1),data',wl)'; % Resampling to the actual spectral sampling interval of the measurement

% Thetas = input('Illumination zenith angle Theta_s (in degree) ?');

% P0=[  N   Cab Car Cbrown  Cw      Cm      Thetai  Bspec]
P0=[    1.5 50  10  0       0.01    0.01    20      0];
LB=[    1   0   0   0       0.0005  0.001   0       -0.2];
UB=[    3   100 30  5       0.1     0.02    90      0.6];

% The "TolFun" parameter has to be tuned depending on signal SNR (low
% tolerance for low SNRs and vice versa)
options = optimset('Algorithm','trust-region-reflective','TolFun',1e-4); 
sol= lsqcurvefit(@procosine,P0,wl,rmes,LB,UB,options,Thetas,data);

rsim=procosine(sol,wl,Thetas,data);

figure,p=plot(wl,rmes,wl,rsim);
set(p(1),'LineWidth',2.5,'Color',[0.75,0.75,0.99])
set(p(2),'LineWidth',1.5,'Color',[0,0,0.7],'LineStyle',':')
xlabel('Wavelength (nm)','FontSize',11,'FontWeight','Bold')
ylabel('Pseudo Bidirectional Reflectance Factor','FontSize',12,'FontWeight','Bold')
title(['N = ',num2str(sol(1),'%4.2f'),'  C_{ab} = ',num2str(sol(2),'%4.1f'),...
    '  C_{cx} = ',num2str(sol(3),'%4.2f'),'  C_{bp} = ',num2str(sol(4),'%4.1f'),...
    '  C_{w} = ',num2str(sol(5),'%7.5f'),'  C_{dm} = ',num2str(sol(6),'%7.5f'),...
    '  \theta_{i} = ',num2str(sol(7),'%4.1f'),'  b_{spec} = ',num2str(sol(8),'%4.2f')],'FontWeight','Bold')
axis([min(wl) max(wl) 0 1])
legend('Measured Pseudo-BRF','Inverted Pseudo-BRF','Location','SouthOutside','Orientation','Horizontal')
