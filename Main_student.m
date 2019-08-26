% Instrumentation and multimodality imaging
% Homework 2
% 2015/2016
% Main (student version)
% 
close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      1. X-Ray simulation          %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_tissues=2;

%%------ fill in with the properties of the mouse data
pixel_size=0.04 ;  % cm 
data_size=[380 380 190] ;

%%------ fill in with the densities of bone and soft tissue at http://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
density_soft=1;
density_bone=1.92;


%% Reading density maps
density_map=zeros(data_size(1),data_size(2),data_size(3),n_tissues);
fid = fopen('digimouse_soft_tissue_sq.img');
density_map(:,:,:,1) = reshape(fread(fid, data_size(1)*data_size(2)*data_size(3),'float'),data_size);
fclose(fid);

fid = fopen('digimouse_bone_sq.img');
density_map(:,:,:,2) = reshape(fread(fid, data_size(1)*data_size(2)*data_size(3),'float'),data_size);
fclose(fid);

figure(1)
subplot(n_tissues,1,1);
imagesc(density_map(:,:,data_size(3)/2,1));
colormap(gray);
axis image
colorbar
subplot(n_tissues,1,2);
imagesc(density_map(:,:,data_size(3)/2,2));
colormap(gray);
axis image
colorbar

%%%% Reading spectra
I0=10^6;
[energy,polySpectrum]=textread('spec_norm_105kVp_1mmAl.txt','%f %f');
polySpectrum105=I0*(polySpectrum/sum(polySpectrum));

[energy,polySpectrum]=textread('spec_norm_65kVp_1mmAl.txt','%f %f');
polySpectrum65=I0*(polySpectrum/sum(polySpectrum));

n_energy=size(energy,1);

%%%% Reading mass attenuation coefficients
MAC=zeros(n_energy,n_tissues);
aux=textread('mac_water.txt', '%f');
MAC(:,1) = squeeze(aux);
aux=textread('mac_bone.txt', '%f');
MAC(:,2) = squeeze(aux);
   
proj_angle=90;
[Raw_Xray_90, Attenuation_image_90] = X_ray_simulator_student(data_size, polySpectrum105, proj_angle, n_tissues, MAC, density_map, pixel_size);

figure(2);
subplot(2,1,1)
imagesc(imrotate(Raw_Xray_90,-90));
colorbar
subplot(2,1,2)
imagesc(imrotate(Attenuation_image_90,-90));
axis image;
colorbar
colormap(gray);

proj_angle=0;
[Raw_Xray_0, Attenuation_image_0] = X_ray_simulator_student(data_size, polySpectrum105, proj_angle, n_tissues, MAC, density_map, pixel_size);

figure(3);
subplot(2,1,1)
imagesc(imrotate(Raw_Xray_0,-90));
colorbar
subplot(2,1,2)
imagesc(imrotate(Attenuation_image_0,-90));
axis image;
colorbar
colormap(gray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%           2. Dual energy calibration       %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating a calibration phantom of the same size of the mouse data set composed of two materials
air=50;
[cal_phantom,mask_soft,mask_bone]=GeneratePhantomSq_student(data_size,air,[density_soft, density_bone]);

figure(4);
imagesc(cal_phantom(:,:,data_size(3)/2));
axis image;
colorbar
colormap(gray);

%------------ Prepare input variable density_map_phantom
density_map_phantom=zeros(data_size(1),data_size(2),data_size(3),n_tissues);

for i=1:data_size(3)
    density_map_phantom(:,:,i,1) = mask_soft;
    density_map_phantom(:,:,i,2) = mask_bone;
end
% Obtaining low energy and high energy projections of the calibration phantom (simulating acquisition) 
proj_angle=0;

[Raw_XrayLow, PhantomLowEnergyProj] = X_ray_simulator_student(data_size, polySpectrum65, proj_angle, n_tissues, MAC, density_map_phantom, pixel_size);

[Raw_XrayHigh, PhantomHighEnergyProj] = X_ray_simulator_student(data_size, polySpectrum105, proj_angle, n_tissues, MAC, density_map_phantom, pixel_size);

figure(5);
subplot(2,1,1)
imagesc(imrotate(PhantomLowEnergyProj,-90));
colorbar
subplot(2,1,2)
imagesc(imrotate(PhantomHighEnergyProj,-90));
axis image;
colorbar
colormap(gray);

% Calculating calibration coefficients
[coeff_soft coeff_bone]=DualEnergyCalibration_student(PhantomLowEnergyProj,PhantomHighEnergyProj, cal_phantom, [density_soft, density_bone],pixel_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      3. Dual Energy Subtraction   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Obtaning low energy and high energy projections of the mouse (simulating acquisition)
proj_angle=90;

[Raw_Xray, MouseLowEnergyProj] = X_ray_simulator_student(data_size, polySpectrum65, proj_angle, n_tissues, MAC, density_map, pixel_size);

[Raw_Xray, MouseHighEnergyProj] = X_ray_simulator_student(data_size, polySpectrum105, proj_angle, n_tissues, MAC, density_map, pixel_size);

% Performing dual energy subtraction
[soft_tissue,bone]=DualEnergySubtraction_student(MouseLowEnergyProj,MouseHighEnergyProj, coeff_soft,coeff_bone);

figure(6);
subplot(2,1,1)
imagesc(imrotate(bone,-90))
colorbar
subplot(2,1,2)
imagesc(imrotate(soft_tissue,-90))
axis image;
colorbar
colormap(gray);