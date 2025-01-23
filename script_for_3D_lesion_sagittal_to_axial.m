%% MRI Lesion Analysis Script
% This script processes MRI lesion data by coregistering, analyzing, and
% calculating projections and affected regions.
% 1. Loading lesion and spinal cord segmentation files.
% 2. Coregistering lesions to the spinal cord template PAM50.
% 3. Projecting lesions and spinal cord segmentations onto axial, sagittal, and coronal planes.
% 4. Calculating the percentage of affected regions in each plane.
% 5. Analyzing masks to determine affected areas within defined regions.
% 6. Exporting results to an Excel file for further analysis.


%% Initialize SPM 
clc;
clear;

% Specify the parent directory containing the MRI data
parentDirectory = ''; 
cd(parentDirectory); 

existingData = {'folderName','Session','axial','sagittal','coronal','MaskPath'};


%Directory containing spinal cord masks for lesion overlap analysis
masksDirectory = '';  

% Get subdirectories within the parent directory
files = dir(parentDirectory);
subdirectories = files([files.isdir] & ~ismember({files.name}, {'.', '..'}));

%% Iterate through each subdirectory
for k = 1:2%length(subdirectories)
    cd(parentDirectory);
    folderName = subdirectories(k).name;
    cd(folderName);
    currentDirectory = pwd;
   
    % Load Lesion and Spinal Cord (SC) segmentation data
    filesInSubdir = dir(currentDirectory);
    scans = filesInSubdir([filesInSubdir.isdir] & ~ismember({filesInSubdir.name}, {'.', '..'}));

    for j = 1:length(scans)
        cd(currentDirectory);  
        scanName = scans(j).name;
        cd(scanName);
            
        % Obtain paths to the lesion and spinal cord segmentation files in
        % PAM 50 space
        lesion = spm_select('FPListRec', pwd, '^volumesoft.*_norm.nii');
        lesionPath = spm_select('FPListRec', pwd, '^volumesoft.*_norm.nii');
        scSegPath = spm_select('FPListRec', pwd, '.*pred-sc_norm.nii');
        
       if isempty(lesionPath)
            disp('Finished');
            continue; % Move to the next iteration of the outer loop
       end
        
        disp(folderName);
        disp(lesionPath);

        %% Interpolate the lesion and volume to match the resolution of the template (0.5 0.5 0.5 mm)
        % %Coregister lesion and spinal cord segmentation
        matlabbatch{1}.spm.spatial.coreg.write.ref = {[scSegPath ',1']};
        matlabbatch{1}.spm.spatial.coreg.write.source = {[lesionPath ',1']};
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;%trilinear interprolation
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

        spm_jobman('run', matlabbatch);
        clear matlabbatch;

        % Read the coregistered lesion volume and spinal cord segmentation
        [lesionPath,file_name]=fileparts(lesionPath);
        lesionVolume = spm_vol(fullfile(lesionPath,['r' file_name '.nii']));
        lesionVolumeMatrix = spm_read_vols(lesionVolume);
        lesionVolumeMatrix(isnan(lesionVolumeMatrix)) = 0;
        lesionVolumeMatrix(lesionVolumeMatrix<0.16) = 0;
        lesionVolumeMatrix(lesionVolumeMatrix>1) = 1;
        lesionBinary = lesionVolumeMatrix; %> 0.16; % Perform trilinear interpolation. Theoretical minimum value is 1/6 => 1 voxel 1 rest 0.

        scSegmentationVolume = spm_vol(scSegPath);
        [scSegPath,file_name_2]=fileparts(scSegPath);
        scSegMatrix = spm_read_vols(scSegmentationVolume);
        scSegMatrix(isnan(scSegMatrix)) = 0;
        scSegMatrix(scSegMatrix<0.16) = 0;
        scSegMatrix(scSegMatrix>1) = 1;
        scSegBinary = scSegMatrix; %> 0.16; % Perform trilinear interpolation. Theoretical minimum value is 1/6 => 1 voxel 1 rest 0.

        %% Apply restriction: lesion has to be within the spinal cord
        lesionBinaryRestricted = lesionBinary .* scSegBinary;

        %% Perform shadowing to create projections
        % Calculate voxel counts and indicators for slices containing lesions
        [V1, V2] = findLesionSlices(lesionBinaryRestricted);

        % Project lesion and spinal cord onto axial, sagittal, and coronal planes
        [axialProjection, sagittalProjection, coronalProjection] = projectLesionOntoPlanes(lesionBinaryRestricted, V1, V2);
        [axialProjectionCord, sagittalProjectionCord, coronalProjectionCord] = projectLesionOntoPlanes(scSegBinary, V1, V2);

        %% Save projection planes
        
        x = size(lesionBinaryRestricted,1);
        y = size(lesionBinaryRestricted,2);
        ii=V1:V2;
        z=size((ii-V1)+1,2);
        
        saveProjectionPlanes(lesion,x,y,z,axialProjection, sagittalProjection, coronalProjection, axialProjectionCord, sagittalProjectionCord, coronalProjectionCord, folderName, lesionPath);

        %% Calculate percentage affected and output results
        [per100Axial, per100Sagittal, per100Coronal] = calculatePercentageAffected(axialProjection, sagittalProjection, coronalProjection, axialProjectionCord, sagittalProjectionCord, coronalProjectionCord);

        output1 = {lesionPath, 'axial', per100Axial};
        output2 = {lesionPath, 'sagittal', per100Sagittal};
        output3 = {lesionPath, 'coronal', per100Coronal};


        disp(output1);
        disp(output2);
        disp(output3);


        % Store data for further analysis or reporting
        newData = {folderName, scanName,num2str(per100Axial), num2str(per100Sagittal), num2str(per100Coronal), 'spinal_cord'};
        existingData = [existingData; newData];

       %% Mask Analysis
        % List all files in the masks directory
        maskFiles = dir(fullfile(masksDirectory, '*.nii'));

    % Iterate through each mask file
        for i = 1:length(maskFiles)
        maskPath = fullfile(masksDirectory, maskFiles(i).name);
        mask = maskFiles(i).name;

        MaskVolume = spm_vol(maskPath);
       [maskPath,file_name_3] = fileparts(maskPath);
        MaskMatrix = spm_read_vols(MaskVolume);
        MaskMatrix(isnan(MaskMatrix)) = 0;
        maskBinary = MaskMatrix; %> 0.16; % Perform trilinear interpolation. Theoretical minimum value is 1/6 => 1 voxel 1 rest 0.
        
        maskBinaryRestricted = lesionBinary .* scSegBinary.* maskBinary;
        TotalmaskBinaryRestricted = scSegBinary.* maskBinary;
       [axialProjectionMask, sagittalProjectionMask, coronalProjectionMask, TotalaxialProjectionMask,TotalsagittalProjectionMask,TotalcoronalProjectionMask] = projectMaskOntoPlanes(maskBinaryRestricted, TotalmaskBinaryRestricted, V1, V2);
       
       %% Save Mask Projections to Nifti
       
        x = size(maskBinaryRestricted,1);
        y = size(maskBinaryRestricted,2);
        ii=V1:V2;
        z=size((ii-V1)+1,2);
        
      saveProjectionMasks(axialProjectionMask, sagittalProjectionMask, coronalProjectionMask, TotalaxialProjectionMask , TotalcoronalProjectionMask ,TotalsagittalProjectionMask, lesionPath,lesion, mask,x ,y,z);
       
           %% Calculate Mask percentage affected and output results
       [per100AxialMask, per100SagittalMask, per100CoronalMask] = calculatePercentageAffectedMask(axialProjectionMask, sagittalProjectionMask, coronalProjectionMask,TotalaxialProjectionMask,TotalsagittalProjectionMask,TotalcoronalProjectionMask);
        output6 = {lesionPath, 'axial', per100AxialMask};
        output7 = {lesionPath, 'sagittal', per100SagittalMask};
        output8 = {lesionPath, 'coronal', per100CoronalMask};

        disp(output6);
        disp(output7);
        disp(output8);


        % Store data for further analysis or reporting
       
        newData = {folderName,scanName, num2str(per100AxialMask),num2str(per100SagittalMask), num2str(per100CoronalMask), mask};
        existingData = [existingData; newData];

        end
    end
end

% Save data to Excel file
outputExcelPath = 'Output.xlsx';
xlswrite(outputExcelPath, existingData);

%% Function to Find Lesion Slices
function [V1, V2] = findLesionSlices(lesionBinaryRestricted)
    % Initialize variables
    indicater = zeros(size(lesionBinaryRestricted,3),1);

    % Identify slices containing lesions
    for i = 1:size(lesionBinaryRestricted,3)
        voxelsLesionOnSlice = sum(sum(lesionBinaryRestricted(:,:,i)));
        if voxelsLesionOnSlice > 0
            indicater(i) = 1;
        end  
    end

    % Find first and last slices containing lesions
    V1 = find(indicater, 1, 'first');
    V2 = find(indicater, 1, 'last');
end

%% Function to Project Lesions onto Planes
function [axialProjection, sagittalProjection, coronalProjection] = projectLesionOntoPlanes(lesionBinaryRestricted, V1, V2)
    % Initialize projection matrices
    axialProjection = zeros(size(lesionBinaryRestricted,1), size(lesionBinaryRestricted,2));
    sagittalProjection = zeros(size(lesionBinaryRestricted,1),((V2-V1)+1));
    coronalProjection = zeros(size(lesionBinaryRestricted,2),((V2-V1)+1));

    % Perform shadowing for axial projection
    for x = 1:size(lesionBinaryRestricted,1)
        for y = 1:size(lesionBinaryRestricted,2)
            axialProjection(x,y) = sum(lesionBinaryRestricted(x,y,V1:V2));
        end
    end

     % Perform shadowing for coronal and sagittal projections
    for ii=V1:V2
        z=(ii-V1)+1;
        %Perform shadowing for coronal
        for x=1:size(lesionBinaryRestricted,1)
            coronalProjection(x,z)=sum(lesionBinaryRestricted(x,:,ii));
        end

        %Perform shadowing for sagittal
        for y=1:size(lesionBinaryRestricted,2)
            sagittalProjection(y,z)=sum(lesionBinaryRestricted(:,y,ii));
        end
    end

   
end

%% Function to Save Projection Planes
function saveProjectionPlanes(lesion,x,y,z,axialProjection, sagittalProjection, coronalProjection, axialProjectionCord, sagittalProjectionCord, coronalProjectionCord, folderName, lesionPath)
    % Create a new directory for saving projection planes
    saveDirectory = fullfile((lesionPath), 'Projection_Planes');
    if ~exist(saveDirectory, 'dir')
        mkdir(saveDirectory);
    end

    % Save axial projection
    
        temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, y, 1];
        HDR_info.fname = ([saveDirectory '\' folderName '_axial_projection.nii']);
        spm_write_vol(HDR_info, axialProjection);

    % Save sagittal projection
        temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [y, z, 1];
        HDR_info.fname = ([saveDirectory '\' folderName '_sagittal_projection.nii']);
        spm_write_vol(HDR_info, sagittalProjection);

    % Save coronal projection
        temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, z, 1];
        HDR_info.fname = ([saveDirectory '\' folderName '_coronal_projection.nii']);
        spm_write_vol(HDR_info, coronalProjection);

    % Save axial projection of spinal cord
      temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, y, 1];
        HDR_info.fname = ([saveDirectory '\' folderName '_axial_cord_projection.nii']);
        spm_write_vol(HDR_info, axialProjectionCord);


    % Save sagittal projection of spinal cord
     temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [y, z, 1];
        HDR_info.fname = ([saveDirectory '\' folderName '_sagittal_cord_projection.nii']);
        spm_write_vol(HDR_info, sagittalProjectionCord);

    % Save coronal projection of spinal cord
         temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, z, 1];
        HDR_info.fname = ([saveDirectory '\' folderName '_coronal_cord_projection.nii']);
        spm_write_vol(HDR_info, coronalProjectionCord);
  
end

%% Function to Calculate Percentage Affected
function [per100Axial, per100Sagittal, per100Coronal] = calculatePercentageAffected(axialProjection, sagittalProjection, coronalProjection, axialProjectionCord, sagittalProjectionCord, coronalProjectionCord)
    % Calculate percentage of affected voxels in each projection plane
    totalAxial = sum((axialProjection(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    totalSagittal = sum((sagittalProjection(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    totalCoronal = sum((coronalProjection(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    
    totalAxialCord = sum((axialProjectionCord(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    totalSagittalCord = sum((sagittalProjectionCord(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    totalCoronalCord = sum((coronalProjectionCord(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    
    % Calculate percentages
    per100Axial = (totalAxial / totalAxialCord) * 100;
    per100Sagittal = (totalSagittal / totalSagittalCord) * 100;
    per100Coronal = (totalCoronal / totalCoronalCord) * 100;
    

    
end


%% Function to Analyze Masks
function  [axialProjectionMask, sagittalProjectionMask, coronalProjectionMask, TotalaxialProjectionMask,TotalsagittalProjectionMask,TotalcoronalProjectionMask] = projectMaskOntoPlanes(maskBinaryRestricted, TotalmaskBinaryRestricted, V1, V2);

    axialProjectionMask = zeros(size(maskBinaryRestricted,1), size(maskBinaryRestricted,2));
    sagittalProjectionMask = zeros(size(maskBinaryRestricted,1),((V2-V1)+1));
    coronalProjectionMask = zeros(size(maskBinaryRestricted,2),((V2-V1)+1));

    % Perform shadowing for axial projection
    for x = 1:size(maskBinaryRestricted,1)
        for y = 1:size(maskBinaryRestricted,2)
            axialProjectionMask(x,y) = sum(maskBinaryRestricted(x,y,V1:V2));
        end
    end

     % Perform shadowing for coronal and sagittal projections
    for ii=V1:V2
        z=(ii-V1)+1;
        %Perform shadowing for coronal
        for x=1:size(maskBinaryRestricted,1)
            coronalProjectionMask(x,z)=sum(maskBinaryRestricted(x,:,ii));
        end

        %Perform shadowing for sagittal
        for y=1:size(maskBinaryRestricted,2)
            sagittalProjectionMask(y,z)=sum(maskBinaryRestricted(:,y,ii));
        end
    end

    TotalaxialProjectionMask = zeros(size(TotalmaskBinaryRestricted,1), size(TotalmaskBinaryRestricted,2));
    TotalsagittalProjectionMask = zeros(size(TotalmaskBinaryRestricted,1),((V2-V1)+1));
    TotalcoronalProjectionMask = zeros(size(TotalmaskBinaryRestricted,2),((V2-V1)+1));

    % Perform shadowing for axial projection
    for x = 1:size(TotalmaskBinaryRestricted,1)
        for y = 1:size(TotalmaskBinaryRestricted,2)
            TotalaxialProjectionMask(x,y) = sum(TotalmaskBinaryRestricted(x,y,V1:V2));
        end
    end

     % Perform shadowing for coronal and sagittal projections
    for ii=V1:V2
        z=(ii-V1)+1;
        %Perform shadowing for coronal
        for x=1:size(TotalmaskBinaryRestricted,1)
            TotalcoronalProjectionMask(x,z)=sum(TotalmaskBinaryRestricted(x,:,ii));
        end

        %Perform shadowing for sagittal
        for y=1:size(TotalmaskBinaryRestricted,2)
            TotalsagittalProjectionMask(y,z)=sum(TotalmaskBinaryRestricted(:,y,ii));
        end
    end
   
    end

%% Function to save Mask Files
function saveProjectionMasks(axialProjectionMask, sagittalProjectionMask, coronalProjectionMask, TotalaxialProjectionMask , TotalcoronalProjectionMask ,TotalsagittalProjectionMask, lesionPath,lesion, mask, x, y,z)

% Create a new directory for saving projection planes
    saveDirectory = fullfile((lesionPath), 'Projection_Masks');
    if ~exist(saveDirectory, 'dir')
        mkdir(saveDirectory);
    end
    
    % Save axial projection
        temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, y, 1];       
        HDR_info.fname = ([saveDirectory '\' mask '_axial_projection.nii']);
        spm_write_vol(HDR_info, axialProjectionMask);

        
    % Save sagittal projection
       temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [y, z, 1];        
        HDR_info.fname = ([saveDirectory '\' mask '_sagittal_projection.nii']);
        spm_write_vol(HDR_info, sagittalProjectionMask);
       

    % Save coronal projection
       temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, z, 1];        
        HDR_info.fname = ([saveDirectory '\' mask '_coronal_projection.nii']);
        spm_write_vol(HDR_info, coronalProjectionMask);
      
    
     % Save axial projection
       temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, y, 1];
        %[file_path,file_name_m] = fileparts(maskPath);
        HDR_info.fname = ([saveDirectory '\' mask '_total_axial_projection.nii']);
        spm_write_vol(HDR_info, TotalaxialProjectionMask);
      
     % Save sagittal projection
        temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [y, z, 1];       
        HDR_info.fname = ([saveDirectory '\' mask '_total_sagittal_projection.nii']);
        spm_write_vol(HDR_info, TotalsagittalProjectionMask);
      
     % Save coronal projection
      temp1=spm_vol(lesion);
        HDR_info = temp1;
        HDR_info.dim = [x, z, 1];
        HDR_info.fname = ([saveDirectory '\' mask '_total_coronal_projection.nii']);
        spm_write_vol(HDR_info, TotalcoronalProjectionMask);
        


end

    %% Function to Calculate Percentage Affected
    
function [per100AxialMask, per100SagittalMask, per100CoronalMask] = calculatePercentageAffectedMask(axialProjectionMask, sagittalProjectionMask, coronalProjectionMask,TotalaxialProjectionMask,TotalsagittalProjectionMask,TotalcoronalProjectionMask)
    % Calculate percentage of affected voxels in each projection plane
    AxialMask = sum((axialProjectionMask(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    SagittalMask = sum((sagittalProjectionMask(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    CoronalMask = sum((coronalProjectionMask(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    
    
    totalAxialMask = sum((TotalaxialProjectionMask(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    totalSagittalMask = sum((TotalsagittalProjectionMask(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    totalCoronalMask = sum((TotalcoronalProjectionMask(:)>0)); %Get rid of the information about how many voxels are present in the third dimension. 
    
    % Calculate percentages
    per100AxialMask = (AxialMask / totalAxialMask) * 100;
    per100SagittalMask = (SagittalMask / totalSagittalMask) * 100;
    per100CoronalMask = (CoronalMask / totalCoronalMask) * 100;
    
end
