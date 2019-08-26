% Instrumentation and multimodality imaging
% Homework 1
% Gabriel Rodríguez Maroto && Andrés Martínez Mora

clear all
close all
clc

%% Mass attenuation coefficient of water (similar to soft tissue).
%------------- CHANGE THE PATH WITH YOUR DIRECTORY---------------
path_mac='.\mac_water.txt';
[Energy,MAC] = textread(path_mac, '%f %f');  % auxiliary vector for reading the attenuation coefficient
KVp=110;
MAC=MAC(1:round(KVp));
Energy=Energy(1:round(KVp));
density=1; %density of water(gr/cm3)
I0=10^6; %number of photons (intensity)
Thickness=3; %cm
      
%% Polychromatic spectrum
%------------- CHANGE THE PATH WITH YOUR DIRECTORY---------------
path_polySpectra='.\Spectra_110KVp.txt'; 
[energy,polySpectra]=textread(path_polySpectra,'%f %f');
energy=energy(1:round(KVp));
polySpectra=polySpectra(1:round(KVp));
polySpectra=I0*polySpectra;
    
%----CALCULATE THE SPECTRUM AFTER TRAVERSING 3 CM OF WATER: IMPLEMENT THE BEER-LAMBERT LAW FOR POLYCHROMATIC SPECTRA USING THE DENSITY AND THE MASS ATTENUATION COEFFCIENT OF WATER ---
mu_water=density*MAC;
polySpectraFilt=polySpectra.*exp(-mu_water*Thickness);
     
figure
plot(energy,polySpectra)
title('Polyenergetic Spectrum')
xlabel('Energy (keV)')
ylabel('Number of Photons (Intensity)')
hold on
plot(energy,polySpectraFilt,'r')
legend('No filtered','Filtered')
 
%% Equivalent monoenergetic spectrum
%--------- WRITE THE MEAN ENERGY OF THE POLYCHROMATIC SPECTRUM EXTRACTED FROM THE SPEKTR INTERFACE -------------
MeanEnergy=sum(dot(Energy,polySpectra))/sum(polySpectra);  %Mean energy of the polyenergetic spectrum
    
MAC_MeanEnergy=MAC(floor(MeanEnergy));%Mass attenuation coefficient of the soft tissue at the mean energy    
             
monoSpectra=[zeros(1,floor(MeanEnergy)-1) I0 zeros(1,energy(end)-floor(MeanEnergy))];  %Monoenergetic spectrum formed by zeros and I0 in the position of MeanEnergy
    

%-------- IMPLEMEMT THE LAMBERT-BEER LAW FOR THE MONOENERGETIC SPECTRUM (ONE VALUE) -----------
monoValueFilt=monoSpectra(floor(MeanEnergy))*exp(-MAC_MeanEnergy*density*Thickness); % Use only the attenuation coefficient corresponding to the energy value of the monochromatic spectrum
   
monoSpectraFilt=[zeros(1,floor(MeanEnergy)-1) monoValueFilt zeros(1,energy(end)-floor(MeanEnergy))];


        
figure
plot(energy,monoSpectra)
title('Monoenergetic Spectrum')
xlabel('Energy (keV)')
ylabel('Number of Photons (Intensity)')
hold on
plot(energy,monoSpectraFilt,'r')
legend('No filtered','Filtered')

%---- FIND THE HVL (ALUMINIUM) OF THE POLYCHROMATIC SPECTRA WITH SIMULATIONS USING THE BEER LAW -----
% To find the HVL of a polychromatic radiation using simulations, we have
% to measure the number of photons trespassing a plate with a certain
% thickness in comparison to the number of photons detected without plate.
% For the thickness where the photons have been reduced to a half, that
% will be the HVL.

% Obtain the linear attenuation coefficient of aluminium from SPEKTR
path_aluminium='.\mu_aluminium.txt';
[energy, mu_aluminium]=textread(path_aluminium,'%f %f');
mu_aluminium=mu_aluminium(1:KVp);
% We will try with thicknesses going from 0 to 10 cm every 0.05 cm
test_thickness=0:0.00001:1;
cont=1; % Counter to go through all the indexes for the thickness when applying the Beer-Lambert's law


% The total number of photons received per thickness will be the sum of all
% the photons with different energies. As we are working with different
% energies between 0 and 110 kVp, the attenuation coefficient for each
% thickness will also be the sum of all the attenuation coefficients. 
% We take as incoming number of photons the number of photons provided by
% the polychromatic spectrum.
% Obtain the ratio per thickness
ratio=zeros(1,length(test_thickness));
I=zeros(1,length(test_thickness));
for i=test_thickness
   I(1,cont)=sum(polySpectra.*exp(-mu_aluminium*i));
   ratio(1,cont)=I(1,cont)/sum(polySpectra);
   cont=cont+1;
end
figure
plot(test_thickness*10,ratio), xlabel('Thickness(mm)'), ylabel('Attenuation ratio (I/I0)'), title('Attenuation ratio vs thickness')

% Now the thickness where the attenuation ratio is the
% closest possible to 0.5 will be the HVL of the polychromatic spectrum.
ind=find(abs(ratio-0.5)==min(abs(ratio-0.5)));
HVL_final=test_thickness(ind);
fprintf('The HVL obtained with simulations of the Beer-Lambert law is %.4f mm\n',HVL_final*10);