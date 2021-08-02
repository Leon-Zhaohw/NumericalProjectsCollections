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

function R_profile = PBD_profile( sigma_a,sigma_s,g,eta,r_array)
%sigma_a: absorption coefficient
%sigma_s: scattering coefficient
%g is the mean cosine, will be used to calculate the reduces parameters with similarity theory
%eta: relative index of refraction
%r_array: distances at which the profile is evaluated

%Init reduced parameters, needed for parameters in the blending range
sigmap_s = sigma_s*(1-g);
sigmap_t = sigma_a + sigmap_s;

%Set number of samples for the equiangularly sampled, near range and
%exponentiallly sampled far range, change for different accuracies
num_samples_equi = 5; 
num_samples_exp = 5; 

%Sets the MIS blending area, which is a very small range around sigma_t 
blending_range_min = 0.9*sigmap_t;
blending_range_max = 1.1*sigmap_t;


function [S_phi, S_E] = eval_sample(in_r,t)
%in_r: sample distance   
%t: distance along extended source 

    %in the orthogonal case, the depth of the source is the distance along the source 
    zr = t;

    %Distance to positive and negative source
    dr = sqrt(in_r.^2 + t.^2); 
    dv = sqrt(in_r.^2 + (t+2.*z_b).^2);
    
    kappa = 1-exp(-2*sigmap_t.*(dr+t)); %Equation (12) in [Habel13pbd]
    
    %Extended source Q(t)
    Q = alpha_p.*sigmap_t.*exp(-sigmap_t.*t); % Equation (23) in [Habel13cid]
 
    R_phi = alpha_p./(4.*pi)./D_g.*(exp(-sigma_tr.*dr)./dr - exp(-sigma_tr.*dv)./dv); % Equation (5) in [Habel13pbd] ((26) in [Habel13cid]) without C_phi
    R_E = alpha_p./(4.*pi).*(zr.*(1+sigma_tr.*dr).*exp(-sigma_tr.*dr)./dr.^3 + (zr+2.*z_b).*(1+sigma_tr.*dv).*exp(-sigma_tr.*dv)./dv.^3); % Equation (5) in [Habel13pbd] ((27) in [Habel13cid])without C_E
    
    S_phi = R_phi.*Q.*kappa;
    S_E = R_E.*Q.*kappa;
     
end

%Grosjean diffusion coefficient
D_g = (2.*sigma_a+sigmap_s)./( 3.*(sigma_a+sigmap_s).^2); %Equation (16) in [Habel13cid]

sigma_tr = sqrt(sigma_a./D_g); %Section 4.1 in [Habel13cid]
alpha_p = sigmap_s/sigmap_t;

A_boundary = (1+QC2x3(eta))/(1-QC1x2(eta)); %Section "4.1: Boundary Conditions" in [Habel13cid]
z_b = 2*A_boundary*D_g;

C_phi = 1/4.*(1-QC1x2(eta)); %Table (2) in [Habel13cid]
C_E = 1/2.*(1-QC2x3(eta)); %Table (2) in [Habel13cid]

rn = 0.5;

r_size = numel(r_array);
R_profile = zeros(r_size,1);
for r_count = 1:r_size
% For clarity reasons, we implement all loops explicitly, the following can
% be formulated as operating on a the whole array of r as it is usually
% done in Matlab. 

    rs = r_array(r_count);
    
    %The following code implements equation (11) in [Habel13pbd].
    %This implementation allows to set different number of samples for
    %the near and far range.

    % Linear weight of the balance heuristic to blend between equiangular and exponential sampling
    % if a weight drops to 0, the corresponding samples are not generated 
    weight = linearstep(blending_range_min,blending_range_max,rs); 

    KP_phi_equi = 0; KP_E_equi = 0;
    if rs < blending_range_max  %only calculate equiangular sampling if in near field or in blending range, otherwise contribution is 0
        for j_equi = 0:(num_samples_equi-1)

            x = (j_equi+rn)/num_samples_equi;
            [t_equi, pdf_t_equi] = equiangular_sampling(x, 0, Inf, rs);
            [S_phi_equi, S_E_equi] = eval_sample(rs,t_equi);
            %get exponential pdf of the generated equiangular sample 
            pdf_t_exp  = exponential_sampling_pdfeval(t_equi,sigmap_t);

            %calculate the balance heuristic weighted by the linear weight from the blending area
            w_t_equi_MIS = (1-weight).*num_samples_equi.*pdf_t_equi./((1-weight).*num_samples_equi.*pdf_t_equi +weight.*num_samples_exp.*pdf_t_exp);

            KP_phi_equi = KP_phi_equi + S_phi_equi.*w_t_equi_MIS./pdf_t_equi;
            KP_E_equi = KP_E_equi + S_E_equi.*w_t_equi_MIS./pdf_t_equi;
        end

        KP_phi_equi = KP_phi_equi./num_samples_equi;
        KP_E_equi = KP_E_equi./num_samples_equi;
    end

    KP_phi_exp = 0; KP_E_exp = 0;
    if rs > blending_range_min %only calculate exponential sampling if in far field or blending range, otherwise contribution is 0
        for j_exp = 0:(num_samples_exp-1)

            x = (j_exp+rn)/num_samples_exp;
            [t_exp, pdf_t_exp] = exponential_sampling(x,sigmap_t);
            [S_phi_exp, S_E_exp] = eval_sample(rs,t_exp);
            %get equiangular pdf of the generated exponential sample 
            pdf_t_equi = equiangular_sampling_pdfeval(t_exp, 0, Inf, rs); 

            %calculate the balance heuristic plus the linear weight from the blending area
            w_t_exp_MIS = weight.*num_samples_exp.*pdf_t_exp./((1-weight).*num_samples_equi.*pdf_t_equi +weight.*num_samples_exp.*pdf_t_exp);   

            KP_phi_exp = KP_phi_exp + S_phi_exp.*w_t_exp_MIS./pdf_t_exp;
            KP_E_exp = KP_E_exp + S_E_exp.*w_t_exp_MIS./pdf_t_exp;
        end

        KP_phi_exp = KP_phi_exp./num_samples_exp;
        KP_E_exp = KP_E_exp./num_samples_exp;
    end

    R_profile(r_count) = C_phi.*(KP_phi_equi+KP_phi_exp) + C_E.*(KP_E_equi+KP_E_exp);
       
end
end

%Equiangular sampling only:  

%         KP_phi = 0; KP_E = 0;
%         for j = 0:(num_samples_equi-1)
% 
%             xi = (j+rn)/num_samples_equi;
%             [t_equi, pdf_t_equi] = equiangular_sampling(xi, 0, Inf, rs);    
%             [S_phi_equi, S_E_equi] = EvalSample(rs,t_equi);
% 
%             KP_phi = KP_phi + S_phi_equi./pdf_t_equi;
%             KP_E = KP_E + S_E_equi./pdf_t_equi;
%         end
% 
%         KP_phi = KP_phi./num_samples_equi;
%         KP_E = KP_E./num_samples_equi;
%     
%         R_profile(r_count) = C_phi.*(KP_phi) + C_E.*(KP_E);
     
%Exponential sampling only:
    
%         KP_phi = 0; KP_E = 0;
%         for j = 0:(num_samples_exp-1)
% 
%             xi = (j+rn)/num_samples_exp;
%             [t_exp, pdf_t_exp] = exponential_sampling(xi,sigma_tp);    
%             [S_phi_exp, S_E_exp] = EvalSample(rs,t_exp);
% 
%             KP_phi = KP_phi + S_phi_exp./pdf_t_exp;
%             KP_E = KP_E + S_E_exp./pdf_t_exp;
%         end
% 
%         KP_phi = KP_phi./num_samples_exp;
%         KP_E = KP_E./num_samples_exp;
%         
%         R_profile(r_count) = C_phi.*(KP_phi) + C_E.*(KP_E);
%     else


