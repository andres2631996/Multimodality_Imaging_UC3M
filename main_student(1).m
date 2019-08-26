% Intrumentation and multimodal imaging
% HOMEWORK 3:PET: Attenuation correction using CT images
% Main

close all
clear all
clc

    % Reading PET files in DICOM format  
    pathPET = '.\Data\dicom_20060905_104748083_NAC'; %path to folder with files
    filesPET = dir(pathPET);
    filesPET(1:2) = [];

	% Loading dicom files ordered by image number 
    for i = 1:size(filesPET,1)
        image_name = filesPET(i).name;
        file = [pathPET,'\',image_name];
        dicom_info = dicominfo (file);          %read dicom info to get the order of the images
        z = dicom_info.InstanceNumber;
        Rescale= double([dicom_info.RescaleSlope]);
        PET_volume(:,:,z) = double(dicomread(file)).*Rescale;    %read dicom images       
        clear image_name z
    end
    
%%%%%%%%%%%% Find information within the PET study (common for all slices)
    weight        =  dicom_info.PatientWeight*1000  ;%(g)
    injected_dose =   dicom_info.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;%(bq) 
    pix_size      =  dicom_info.PixelSpacing/10  ;%cm 
    Slice_thick   = dicom_info.SliceThickness/10   ;%cm 
    PET_size      = [dicom_info.Rows dicom_info.Columns dicom_info.NumberOfSlices] ;
    PETunits      =  dicom_info.Units; 
    PET_VoxelSize = [pix_size(1) pix_size(2) Slice_thick]; 
    
    % Visualizing the PET volume
    %figure, imshow3D(PET_volume) 
    figure, imshow(PET_volume(:,:,138),[]);
    figure, imshow(squeeze(PET_volume(PET_size(1)/2,:,:))',[])

    
    % Reading the CT for PET attenuation correction
%%%%%%%%%%%%% Find matrix and voxel size for CTAC volume
% Read the header with the Matlab command textscan to read line by line
% organizing each line by reading a string and three numeric variables that
% are saved into the same vector depending on the field (Voxel Size or
% Matrix Size)
    path_CT='.\Data\CT_files';
    fi=fopen(strcat(path_CT,'\CT_for_ACPET__120KeV_148mA.mhd'));
    mhd_AC=textscan(fi,'%s %f %f %f','Delimiter','=');
    CTAC_size=[];
    CTAC_voxel_size=[];
    for i=2:4
        CTAC_size = [CTAC_size mhd_AC{1,i}(3)] ;
        CTAC_voxel_size= [CTAC_voxel_size mhd_AC{1,i}(2)]  ;
    end
    CTAC_voxel_size=CTAC_voxel_size/10;
    path_CT_AC = '.\Data\CT_files\CT_for_ACPET__120KeV_148mA.raw';
    fileID = fopen(path_CT_AC);
    CT_forAC = reshape(fread(fileID,CTAC_size(1)*CTAC_size(2)*CTAC_size(3),'int16'),CTAC_size);
    fclose(fileID);
    
%%%%%%%% Find voltage for CTAC volume
    voltage       = str2num(path_CT_AC(31:33)) ; %KeV
    
    % Reading the standard CT    
%%%%%%%%%%%%% Find matrix for standard CT volume
% Do the same performed before with the header of the file of CT for
% attenuation correction
    fi=fopen(strcat(path_CT,'\CT_STANDARD_70CM_FOV_3.27MM_120KeV_148mA.mhd'));
    mhd_standard=textscan(fi,'%s %f %f %f','Delimiter','=');
    CT_size=[];
    for i=2:4
        CT_size =[CT_size mhd_standard{1,i}(3)]      ;
    end  
    path_CT_standard = '.\Data\CT_files\CT_STANDARD_70CM_FOV_3.27MM_120KeV_148mA.raw';
    fileID  = fopen(path_CT_standard);
    CT_volume_std = reshape(fread(fileID,CT_size(1)*CT_size(2)*CT_size(3),'int16'),CT_size);
    fclose(fileID); 
    
    % Displaying both CT volumes
    figure
    subplot(1,2,1);
    imshow(CT_volume_std(:,:,(CT_size(3)-1)/2),[]);
    colormap(gray);
    axis image
    colorbar
    subplot(1,2,2);
    imshow(squeeze(CT_forAC(:,:,CTAC_size(3)/2)),[])
    colormap(gray);
    axis image
    colorbar
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  2. Match PET spatial resolution and orientation %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [CT_volume_interp]= PET_size_reformat_student(CT_forAC,CTAC_voxel_size,PET_VoxelSize);

%%%%%%% Check the size and orientation of the new CT volume

    figure, imshowpair(CT_volume_interp(:,:,50),PET_volume(:,:,50),'montage');
    figure, imshowpair(CT_volume_interp(:,:,79),PET_volume(:,:,79),'montage');
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  3.Conversion of HU to attenuation coefficients at 511 keV %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [mu511_map_interp] = CT_energyScaling_student (CT_volume_interp, voltage);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      4. Attenuation correction         %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [PET_corrected] = attenuation_correction_student (PET_volume, mu511_map_interp, PET_size,pix_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%      5. SUV conversion         %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [PET_corrected_SUV] = SUV_conv_student (PET_corrected, injected_dose, weight, PETunits);
    [PET_raw_SUV] = SUV_conv_student (PET_volume, injected_dose, weight, PETunits);
    
    % Displaying the PET volume before and after correction    figure(1)
    subplot(2,2,1);
    imshow(PET_corrected_SUV(:,:,138),[]);
    colormap(gray);
    axis image
    colorbar
    subplot(2,2,2);
    imshow(squeeze(PET_corrected_SUV(PET_size(1)/2,:,:))',[])
    colormap(gray);
    axis image
    colorbar
    
    subplot(2,2,3);
    imshow(PET_raw_SUV(:,:,138),[]);
    colormap(gray);
    axis image
    colorbar
    subplot(2,2,4);
    imshow(squeeze(PET_raw_SUV(PET_size(1)/2,:,:))',[])
    colormap(gray);
    axis image
    colorbar
    
    % Saving both volumes
    path_save = '.\PET_corr_SUV.raw';
    fileID = fopen (path_save,'w');
    fwrite(fileID,PET_corrected_SUV,'float');
    fclose (fileID);
    
    path_save = '.\PET_SUV.raw';
    fileID = fopen (path_save,'w');
    fwrite(fileID,PET_raw_SUV,'float');
    fclose (fileID);