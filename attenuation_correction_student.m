function [PET_corrected] = attenuation_correction_student (PET_volume, mu511_map_interp,PET_size,pix_size)
    
    % Implements PET attenuation correction.
    % Inputs:
    %   PET_volume: original PET volume
    % 	mu511_map_interp: volume with the map of mu values at 511kev interpolated to PET spatial resolution
    % Outputs:
    % 	PET_corrected: PET volume with attenuation correction


    % R = radon(I, theta) returns the Radon transform R of the intensity image I for the angle theta degrees.
    % The Radon transform is the projection of the image intensity along a radial line oriented at a specific angle
    theta = 0:180;
    
%%%%% Compute the attenuation factors 
attenuation_factor=[];
sinogram_pet=[];
PET_size=double(PET_size);
for i=1:size(mu511_map_interp,3)
 attenuation_factor(:,:,i)=exp(-radon(mu511_map_interp(:,:,i)*pix_size(1),theta));
 sinogram_pet(:,:,i)=radon(PET_volume(:,:,i),theta);  
end
%%%%%% Correct the PET volume
PET_corrected=[];
for i=1:size(mu511_map_interp,3)
    PET_corrected(:,:,i)=iradon(sinogram_pet(:,:,i)./attenuation_factor(:,:,i),theta); 
end   
scaleFactor = size(PET_corrected)./PET_size;
[u,v,w] = meshgrid(1:scaleFactor(1):size(PET_corrected,1),1:scaleFactor(2):size(PET_corrected,2),1:scaleFactor(3):size(PET_corrected,3));
PET_corrected=interp3(PET_corrected,u,v,w);
end
