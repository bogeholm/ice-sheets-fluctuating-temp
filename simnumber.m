function [simno] = simnumber(T_want, itm_want, ppfac_want, ...
    T_warming, ppfac, itm_c)
% Convert from (T_want, itm, ppfac) to simulation number in Robinsons
% data

% Index in the vector
[~, IT] = min(abs(T_want - T_warming));
[~, Ipp] = min(abs(ppfac - ppfac_want));
[~, Iit] = min(abs(itm_c - itm_want));
% Simulation number
simno = (IT - 1)*99 + (Ipp - 1)*11 + Iit;

end

