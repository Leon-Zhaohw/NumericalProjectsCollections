%******************************************************************************
%
%Copyright (c) 2013 by Disney-Pixar
%
%Permission is hereby granted to use this software solely for 
%non-commercial applications and purposes including academic or 
%industrial research, evaluation and not-for-profit media
%production.  All other rights are retained by Pixar.  For use 
%for or in connection with commercial applications and
%purposes, including without limitation in or in connection 
%with software products offered for sale or for-profit media
%production, please contact Pixar at tech-licensing@pixar.com.

%******************************************************************************

function R_profile = classical_dipole(sigma_a,sigma_s,g,eta,r_array)
%sigma_a: absorption coefficient
%sigma_s: scattering coefficient
%g is the mean cosine, will be used to calculate the reduces parameters with similarity theory
%eta: relative index of refraction
%r_array: distances at which the profile is evaluated

%calculate reduced scattering parameters
sigmap_s = sigma_s*(1-g);
sigmap_t = sigma_a + sigmap_s;


%classical diffusion coefficient
D = 1 ./ (3 * sigmap_t); %Table (2) in [Habel13cid]

sigma_tr = sqrt(sigma_a./D);
alpha_p = sigmap_s/sigmap_t;

%Boundary conditions from [Jensen01]
F_dr = -1.440 * eta^(-2) + 0.710 * eta^(-1) + 0.668 + 0.0636 * eta;
A_boundary = (1 + F_dr) ./ (1 - F_dr); %Classical boundary, defined in Section 3.2: Boundary Conditions in [Habel13cid]
z_b = 2*A_boundary*D;
zr = 1./sigmap_t;

r_size = numel(r_array);
R_profile = zeros(r_size,1);
for r_count = 1:r_size
% For clarity reasons, we implement all loops explicitly, the following can
% be formulated as operating on a the whole array of r as it is usually
% done in Matlab. 

    rs = r_array(r_count);
    
    dr = sqrt(rs.^2 + zr.^2); %Section 3.2:The Dipole in [Habel13cid] 
    dv = sqrt(rs.^2 + (zr+2.*z_b).^2); %Section 3.2:The Dipole in [Habel13cid] 

    R_E = alpha_p./(4.*pi).*(zr.*(1+sigma_tr.*dr).*exp(-sigma_tr.*dr)./dr.^3 + (zr+2.*z_b).*(1+sigma_tr.*dv).*exp(-sigma_tr.*dv)./dv.^3); %Equation (12) in [Habel13cid]
         
    R_profile(r_count) = R_E; %Equation (20) in [Habel13cid]


end

end

