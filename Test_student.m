
    clear
    close all
    clc

    offset=1024;
    % Number of tranaxial, angular, and axial data samples: IMPLEMENT
    sino_rho = 55;
    sino_phi = 120;      
    det_rows = 28;
   
    
    %Total number of sinograms (direct + crossed): IMPLEMENT
    total_sino= (det_rows)^2 ;

    im_size=[sino_rho,sino_rho];
    % Number of final slices: IMPLEMENT
    % Obtain intermediate slices from adjacent sinograms. In the edge,
    % there is no adjacent sinogram
    n_slices = det_rows*2-1;
    
    % Reconstructed volumes
    rec_image=zeros(sino_rho,sino_rho,n_slices);
    filt_recon_ramp=zeros(sino_rho,sino_rho,n_slices);
    
    %Total rotation angle needed for reconstruction with parallel geometry: IMPLEMENT
    total_angle=180  ;
    
    % read meassured PET data file
    projections_file='./data/CANNABIS_MARZO_2006_Mar01_Acq02_001_001.sin';
    fid = fopen(projections_file, 'rb');
    aux=fseek(fid, offset, 'bof');
    Proj = cast(reshape(fread(fid,sino_rho*sino_phi*total_sino, 'int16', 'ieee-le'),sino_rho,sino_phi,det_rows,det_rows),'double');
    fclose(fid);
    
    % read flood field inverse file
    sensib_file='./data/sensibs.sin';
    fid = fopen(sensib_file, 'rb');
    sensib = cast(reshape(fread(fid,sino_rho*sino_phi*total_sino, 'float', 'ieee-le'),sino_rho,sino_phi,det_rows,det_rows),'double');
    fclose(fid);

    % Sinogram uniformity correction: IMPLEMENT
    sino3d=zeros(sino_rho,sino_phi,det_rows,det_rows);
    % For all sinograms
    sensib=sensib+eps;
    sino3d=Proj./sensib;
   
    % Perform Single Slice Rebinning (SSRB)
    max_diff=7;
    sinogram2D=SSRB_student(sino3d,max_diff, n_slices);
    
    for slice=1:n_slices 
        % The FBP method is applied to the filtered sinogram
        rec_image(:,:,slice)=iradon(sinogram2D(:,:,slice),[1:total_angle/sino_phi:total_angle],'linear','Ram-Lak',im_size(1));
        filt_recon_ramp(:,:,slice)=iradon(sinogram2D(:,:,slice),[1:total_angle/sino_phi:total_angle],'linear','Hann',im_size(1));
    end
    
    figure(1); imagesc(abs(rec_image(:,:,(n_slices-1)/2)));
    colormap(gray)
    figure(2); imagesc(abs(filt_recon_ramp(:,:,(n_slices-1)/2)));
    colormap(gray)
    
    % Saving both volumes
    path_save = '.\RecImage.raw';
    fileID = fopen (path_save,'w');
    fwrite(fileID,rec_image,'float');
    fclose (fileID);
    
    path_save = '.\FiltImage.raw';
    fileID = fopen (path_save,'w');
    fwrite(fileID,filt_recon_ramp,'float');
    fclose (fileID);

    


    
    