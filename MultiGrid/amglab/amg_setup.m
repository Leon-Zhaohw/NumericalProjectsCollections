% AMG_SETUP calls the specified routines for coarse point and intergrid weight selection

%  Ryan McKenzie
%  Department of Computational Sciences
%  University of Kentucky

function amg_setup(level);

amg_globals;

if level ~= COARSEST % if the current level is not the coarsest level
    % call appropriate setup algorithm to select coarse points and get
    % interpolation weights
    
    if SETUP_ALG == RUGE_STUEBEN
        toodeep = RugeStuebenCoarsen(level);
    elseif SETUP_ALG == AMGm
        toodeep = AMGmCoarsen(level);
    elseif SETUP_ALG == SMOOTH_AGGREGATE
        toodeep = SmoothAggregateCoarsen(level);
    elseif SETUP_ALG == AMGe
        toodeep = AMGeCoarsen(level);
    elseif SETUP_ALG == BECK
        toodeep = beck_coarsen(level);
    end
    
    if toodeep==1 % the setup algorithm cannot coarsen to the next level
        COARSEST = level; % make the current level the coarsest
        warning = sprintf('Coarsest level set too deep, resetting to %3d',level)
    end
    
    if SETUP_OPT == AT_ONCE % if we are coursening more than one level
        if level ~= COARSEST % if the current level is not the coarsest level
           amg_setup(level+1); % call setup for next course grid
        end
    end
end
