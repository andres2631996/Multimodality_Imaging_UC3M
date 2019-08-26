% Instrumentation and multimodality imaging
% Homework 2
% 2015/2016
% Dual Energy Calibration (student version)

function [coeff_tissue1 coeff_tissue2]=DualEnergyCalibration(LowEnergyProj,HighEnergyProj, cal_phantom, density_val, pixel_size)

size_data=size(LowEnergyProj);

%% phantom segmentation
phantom_tissue1=cal_phantom*0;
phantom_tissue1(cal_phantom==density_val(1))=cal_phantom(cal_phantom==density_val(1));

phantom_tissue2=cal_phantom*0;
phantom_tissue2(cal_phantom==density_val(2))=cal_phantom(cal_phantom==density_val(2));

%--------- Find the traversed mass thicknesses for each pixel in the projection
mass_thick_tissue1=zeros(1,size_data(2));
mass_thick_tissue2=zeros(1,size_data(2));

mass_thick_tissue1=squeeze(sum(phantom_tissue1,2))*pixel_size; %cm;
mass_thick_tissue2=squeeze(sum(phantom_tissue2,2))*pixel_size; %cm;

mass_thick_tissue1=mass_thick_tissue1(:,size_data(2)/2);
mass_thick_tissue2=mass_thick_tissue2(:,size_data(2)/2);

LowEnergyProj_data=LowEnergyProj(:,size_data(2)/2);
HighEnergyProj_data=HighEnergyProj(:,size_data(2)/2);

% Preparing equation for each pixel value: mass_thick_tissue=CC*M  
Ones=ones(size(LowEnergyProj_data));
M=[Ones,LowEnergyProj_data,HighEnergyProj_data,...
        LowEnergyProj_data.^2,HighEnergyProj_data.^2,...
        LowEnergyProj_data.*HighEnergyProj_data,LowEnergyProj_data.^3, ...
        HighEnergyProj_data.^3];

% Solving the system of equations: mass_thick_tissue=CC*M   
[U,S,V] = svd(M,0);
s=diag(S);
figure, semilogy(s/s(1),'*')
sinv=1./s;
sinv(s/s(1)<1e-6) = 0; % truncated SVD
coeff_tissue1=V*diag(sinv)*transpose(U)*mass_thick_tissue1;
coeff_tissue2=V*diag(sinv)*transpose(U)*mass_thick_tissue2;

end



