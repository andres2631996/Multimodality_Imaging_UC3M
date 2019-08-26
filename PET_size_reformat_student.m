function [CT_volume_interp] = PET_size_reformat_student(CT_volume,CTAC_voxel_size,PET_VoxelSize)

% Interpolation to PET spatial resolution
% Inputs: 
%   CT_volume: CT data with original size
%   CTAC_voxel_size: voxel size of CT volume: [vox_size_x, vox_size_y, vox_size_z]
%   PET_VoxelSize: voxel size of PET volume: [vox_size_x, vox_size_y, vox_size_z]
% Outputs:
%   CT_volume_interp: CT volume interpolated to PET spatial resolution

inputSize=size(CT_volume);

%%%%% Compute scale factor
scaleFactor = PET_VoxelSize./CTAC_voxel_size;

[u,v,w] = meshgrid(1:scaleFactor(1):inputSize(1),1:scaleFactor(2):inputSize(2),1:scaleFactor(3):inputSize(3));

%%%%%%%%%%%%% HERE: CHECK THAT THE FINAL NUMBER OF ELEMENTS IS CORRECT!!!

CT_volume_interp = interp3(CT_volume,u,v,w);

end