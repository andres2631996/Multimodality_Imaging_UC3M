function [PET_SUV] = SUV_conv_student (PET_corrected, injected_dose, weight, PETunits)

    % Converts counts into Standardized Uptake Values (SUV)
    % Inputs:
    % 	PET_volume: PET volume
    % 	Injected_dose: Dose injected to the patient in Bequerels
    % 	Weight: Patient weight in grams
    %   PETunits: Voxel value in PET. If "PROPCPS" then values are in MBq
    % Outputs:
    % 	PET_SUV: PET volume in SUVs
 if strcmp(PETunits,'PROPCPS')==1
     PET_corrected=PET_corrected*10^6;
     PET_SUV=PET_corrected*weight/injected_dose;
 end
    
end