function [q, q_indx] = descendSort(q_dof_val, q_dof)

N = size(q_dof,2);
for l=1:N-1,
    for j =1: N - l,
        if q_dof_val(j) < q_dof_val(j+1)
            temp = q_dof_val(j);
            q_dof_val(j) = q_dof_val(j+1);
            q_dof_val(j+1) = temp;
            temp = q_dof(j);
            q_dof(j) = q_dof(j+1);
            q_dof(j+1) = temp;
        end
    end
end

q = q_dof_val;
q_indx = q_dof;