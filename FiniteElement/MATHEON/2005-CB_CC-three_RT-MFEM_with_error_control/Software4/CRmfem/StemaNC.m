
function M_NC=StemaNC(vertices)
G=[ones(1,3);vertices']\[zeros(1,2);eye(2)];
M_NC=4*det([ones(1,3);vertices'])*G*G';