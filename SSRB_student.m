function sinogram2D=SSRB_student(sino3d, maxring_diff, n_slices)

    tam=size(sino3d);
    sino_rho = tam(1);
    sino_phi  = tam(2);
    det_rows = tam(3);
    
    arr_result = zeros(sino_rho, sino_phi, n_slices);
    arr_contributions = zeros(n_slices,1);
   
    % Loop to go through every 3D sinogram
    for i=1:det_rows     
        for j=1: det_rows
            
            % destination 2D sinogram
            central_z=i+j-1;

            % condition to be added to the same direct sinogram: IMPLEMENT
            if i-j<=maxring_diff
                arr_result(:,:,central_z) = arr_result(:,:,central_z) + sino3d(:,:,i,j);
                arr_contributions(central_z)=arr_contributions(central_z)+1;
            end
        end
    end
    
    % Normalize for the number of sinograms used to create each
    % direct sinogram: IMPLEMENT
    for i=1:n_slices
        sinogram2D(:,:,i)= arr_result(:,:,i)/arr_contributions(i)  ;
    end
end
