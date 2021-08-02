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

function R_profile = single_scattering_profile(sigma_a,sigma_s,g,eta,r_array)
%sigma_a: absorption coefficient
%sigma_s: scattering coefficient
%g is the mean cosine, will be used to calculate the reduces parameters with similarity theory
%eta: relative index of refraction
%r_array: distances at which the profile is evaluated
%set number of samples

num_samples = 5;

sigma_t = sigma_a + sigma_s;
rn = 0.5;

r_size = numel(r_array);
R_profile = zeros(r_size,1);
for r_count = 1:r_size
% For clarity reasons, we implement all loops explicitly, the following can
% be formulated as operating on a the whole array of r as it is usually
% done in Matlab.
 
    rs = r_array(r_count);

    for j = 0:(num_samples-1)

        xi = (j+rn)/num_samples;

        critangle = pi/2 - asin(1/eta);
        if(imag(critangle) ~= 0)
            critangle = 0;
        end
        %start integration from critical angle to infinity
        [t, pdf_t] = equiangular_sampling(xi, rs*tan(critangle), Inf, rs);
        
        % distance to source
        dp = sqrt(rs.^2 + t.^2);

        %not really a Kienle-Pattersen variable, named for consistency with multi scattering implementation
        KP = sigma_s.*HG(g,cos(pi-atan(rs./t))).*exp(-sigma_t.*t).*exp(-sigma_t.* dp)./dp.^2 .*(t./dp).*fresneltrans(cos(pi/2 -atan(t./rs)),1/eta);
        R_profile(r_count) = R_profile(r_count) + KP./pdf_t;

    end

    R_profile(r_count) = R_profile(r_count)./num_samples;


end

end


