% Instrumentation and multimodality imaging
% Homework 2
% 2015/2016
% Generate Square Phantom (student version)

function [cal_phantom,mask_soft,mask_bone]=GeneratePhantomSq_student(matrix_size,air, densities)

cal_phantom=zeros(matrix_size(1),matrix_size(2),matrix_size(3));
phantom_slice=zeros(matrix_size(1),matrix_size(2));
mask_soft=zeros(matrix_size(1),matrix_size(2));
mask_bone=zeros(matrix_size(1),matrix_size(2));

%---------- create a phantom with the shape of a square prism, half
%---------- density(1) and half density (2)

for i=air:matrix_size(1)-air
    for j=air:matrix_size(2)-air
        if i>=j
            phantom_slice(i,j)=densities(2);
            mask_bone(i,j)=densities(2);
        else
            phantom_slice(i,j)=densities(1);
            mask_soft(i,j)=densities(1);
        end
    end
end

for i=1:matrix_size(3)
    cal_phantom(:,:,i)=phantom_slice;
end
end
