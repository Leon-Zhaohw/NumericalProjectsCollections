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

function R_profile = improved_dipole( sigma_a,sigma_s,g,eta,r_array)
%sigma_a: absorption coefficient
%sigma_s: scattering coefficient
%g is the mean cosine, will be used to calculate the reduces parameters with similarity theory
%eta: relative index of refraction
%r_array: distances at which the profile is evaluated

%calculate reduced scatterin parameters
sigmap_s = sigma_s*(1-g);
sigmas_t = sigma_a + sigmap_s;

%Grosjean diffusion coefficient
D_g = (2.*sigma_a+sigmap_s)./( 3.*(sigma_a+sigmap_s).^2); %Equation (16) in [Habel13cid]

sigma_tr = sqrt(sigma_a./D_g);
alpha_p = sigmap_s/sigmas_t;

A_boundary = (1+QC2x3(eta))/(1-QC1x2(eta)); %Section "4.1: Boundary Conditions" in [Habel13cid]
z_b = 2*A_boundary*D_g; %Section 3.2:The Dipole in [Habel13cid] 
zr = 1./sigmas_t; %place source at 1 mean free path

C_phi = 1/4.*(1-QC1x2(eta)); %Table (2) in [Habel13cid]
C_E = 1/2.*(1-QC2x3(eta)); %Table (2) in [Habel13cid]

r_size = numel(r_array);
R_profile = zeros(r_size,1);
for r_count = 1:r_size
% For clarity reasons, we implement all loops explicitly, the following can
% be formulated as operating on a the whole array of r as it is usually
% done in Matlab.  

    rs = r_array(r_count);
    
    dr = sqrt(rs.^2 + zr.^2); %Section 3.2:The Dipole in [Habel13cid] 
    dv = sqrt(rs.^2 + (zr+2.*z_b).^2); %Section 3.2:The Dipole in [Habel13cid] 
    
    %compared to PBD, there is an additional alpha_p to account for the
    %scattering to the point inside the medium which is carried by the
    %extended source in the case of photon beam diffusion
    R_phi = alpha_p^2./(4.*pi)./D_g.*(exp(-sigma_tr.*dr)./dr - exp(-sigma_tr.*dv)./dv); %Equation (21) in [Habel13cid]
    R_E = alpha_p.^2/(4.*pi).*(zr.*(1+sigma_tr.*dr).*exp(-sigma_tr.*dr)./dr.^3 + (zr+2.*z_b).*(1+sigma_tr.*dv).*exp(-sigma_tr.*dv)./dv.^3); %Equation (22) in [Habel13cid]
    
    R_profile(r_count) = C_phi.*(R_phi) + C_E.*(R_E); %Equation (20) in [Habel13cid]

end

end

