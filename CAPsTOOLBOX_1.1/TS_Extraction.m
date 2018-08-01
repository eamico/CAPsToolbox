function TS = TS_Extraction(data,V,seed_mni,seed_radius)

% This function extracts the seed time course from your data.
%  
% usage:
% TS = TS_Extraction(data,V,seed_mni,seed_radius);
%
% Input:
% data - cell array (1,number of subjects). Each cell should contain an 
% NxM matrix for N time points (e.g. functional volumes) in M voxels;
% V - spm vector of structures containing image volume information. Needed 
% as reference for the voxel to world transformation (see spm_vol for info)
% seed_mni- 3D coordinates or array of coordinates (in mm) of your seed 
% centroid(s)  (e.g seed_mni= [0 53 26] or seed_mni= [0 -53 26;0 54 -8;])
% seed_radius = your seed radius (in mm, e.g. seed_radius=6). 

%
% Output: 
% TS - cell array (number of subjects,number of seeds). Each cell contain the seed Time 
% course Mx1 array with M=number of time points (e.g fMRI volumes).
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
% This file is part of CAPsTOOLBOX.  It is Copyright (C) Enrico Amico, 
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
%__________________________________________________________________________
if(iscell(V))
    V=V{1};
end
num_seed=size(seed_mni,1);
TS=cell(length(data),num_seed);
for j=1:num_seed    
    
    for i=1:length(data)
        tmp=data{i};
        fprintf('\n Processing subj %d \n',i);
        [seed_cor] = Comp_Ball2Mask(V.dim,abs(diag(V.mat(1:3,1:3)))',seed_mni(j,:), seed_radius, V);
        seed_ind = sub2ind(V.dim, seed_cor(:,1), seed_cor(:,2), seed_cor(:,3));
        TS{i,j} = mean(tmp(:,seed_ind),2);
    end
    
end