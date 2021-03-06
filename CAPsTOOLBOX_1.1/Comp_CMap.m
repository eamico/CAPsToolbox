function CMap = Comp_CMap(data,V,brind_tot,seed_mni,seed_name,seed_radius)  

% This function calculates seed-voxel correlation maps from your data
%  
%% usage:
%% CMap = Comp_CMap(data,V,brind,seed_mni,seed_name,seed_radius);
%
%% Input:
% data - cell array (1,number of subjects). Each cell should contain an 
% NxM matrix for N time points (e.g. functional volumes) in M voxels;
% V - spm vector of structures containing image volume information. Needed 
% as reference for the voxel to world transformation (see spm_vol for info)
% brind = brain mask, organized as 1xK, with K=number of brain voxels. -Last update EA June 2015: individual brind as input (cell array)           
% Note that brind = find(brainmask(:)), where brainmask is a 3D GrayMatter mask, with 1 in Gray matter voxels and 0 elsewhere 
% seed_mni- 3D coordinates or array of coordinates (in mm) of your seed 
% centroid (e.g seed_mni= [0 -53 26] or seed_mni= [0 -53 26;0 54 -8;])
% seed_name - cell array of string(s) containing the name(s) of the seed(s)
% (e.g. seed_name ={'PCC','MPFC','THA','THA2'} or seed_name ={'PCC'}).  
% seed_radius = your seed radius (in mm, e.g. seed_radius=6). 
%
%% Output: 
% CMap - cell array (number of subjects, number of seeds). Each cell 
% contain the seed-voxel correlation map as a 1xM array with M=number of 
% voxels.It can be saved into an .img/.nii file using the spm_write 
% function (see spm_write help for info). 
%
% Ref: Amico, Enrico, et al. "Posterior Cingulate Cortex-Related 
% Co-Activation Patterns: A Resting State fMRI Study in Propofol-Induced 
% Loss of Consciousness." PloS one 9.6 (2014): e100012.
% 
% Written by EA Sep 2014
% Updated by EA July 2018
%__________________________________________________________________________
%License
%
% This file is part of CAPsTOOLBOX. It is Copyright (C) Enrico Amico. 
% 
% CAPsTOOLBOX is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CAPsTOOLBOX is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.
% __________________________________________________________________________
% 

num_seed=size(seed_mni,1);
CMap=cell(length(data),num_seed);

if(iscell(V))
    V=V{1};
end

for i=1:length(data)
    
    fprintf('\n Processing Subj %d \n',i);
    if(iscell(brind_tot))
        brind=brind_tot{i};
    else
        brind=brind_tot;        
    end
    Y = data{i}(:,brind);
    fprintf('\n Number of NaN in the data %d \n',length(find(isnan(Y))));
    
    for iseed=1:num_seed
        
        [seed_cor] = Comp_Ball2Mask(V.dim,abs(diag(V.mat(1:3,1:3)))',seed_mni(iseed,:), seed_radius, V);
        fprintf('\n Seed %s \n',seed_name{iseed});
        seed_ind = sub2ind(V.dim, seed_cor(:,1), seed_cor(:,2), seed_cor(:,3));
        X = mean(data{i}(:,seed_ind),2);
        fprintf('\n Number of NaN in the seed %d \n',length(find(isnan(X)))); 
        dat = zeros(V.dim);
        tmp = Comp_fastCorr(X, Y); 
        dat(brind) = tmp;
        fprintf('\n Number of NaN in the CMap image %d \n',length(find(isnan(dat))));
        dat(isnan(dat)) = 0;
        CMap{i,iseed} = dat(:)';
    
    end
    
end

return;
