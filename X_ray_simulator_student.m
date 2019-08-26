% Instrumentation and multimodality imaging
% Homework 2
% 2015/2016
% X-Ray Simulator (student version)

function [ Raw_Xray, Attenuation_image ] = X_ray_simulator_student(data_size, polySpectrum, proj_angle, n_tissues, MAC, density_map, pixel_size)

n_energy=size(polySpectrum,1);

% Auxiliar variables
rot_density_map=zeros(data_size(1),data_size(2),data_size(3),n_tissues); % Density map to be rotated
mass_thick=zeros(data_size(1),data_size(3),n_tissues);
mu_energy=zeros(data_size(1),data_size(3),n_energy);

% Resulting images
Raw_Xray=zeros(data_size(1),data_size(3));
Attenuation_image=zeros(data_size(1),data_size(3));


for i=1:n_tissues  
%----- Prepare the projection angle
    if proj_angle~=0
        for j=1:data_size(3) % Rotate in the Z axis to obtain the different projections
            rot_density_map(:,:,j,i)=imrotate(density_map(:,:,j,i),proj_angle);
        end
    else
      rot_density_map=density_map;  
    end
    
%--------- Calculate mass thickness: Integral of density data
        mass_thick(:,:,i)=sum(rot_density_map(:,:,:,i),2)*pixel_size;
 
        
%--------- Update exponent in Beer Lambert Law
    for k=1:n_energy
        mu_energy(:,:,k)=sum(mass_thick(:,:,i)*MAC(k,i),3)+mu_energy(:,:,k);
    end
end

tmp_proj=zeros(data_size(1),data_size(3),n_energy);
for k=1:n_energy
%---------------- Implement Beer Lamber Law    
    tmp_proj(:,:,k)=polySpectrum(k)*exp(-mu_energy(:,:,k));
end

%----------- Obtain resulting images
 Raw_Xray=sum(tmp_proj,3);
 Attenuation_image=-log(Raw_Xray/sum(polySpectrum));

end
