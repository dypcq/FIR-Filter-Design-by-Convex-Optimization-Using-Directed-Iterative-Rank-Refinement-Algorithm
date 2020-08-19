function [maskupper,masklower] = filtermask_v1(freq,publl,pubul,pubv,plbll,plbul,plbv,sll,sul,sv)

% This module produces a mask for filter design.
% 
% freq :frequency vector
% publl:passband upper bound lower (frequency)limit
% pubul:passband upper bound upper (frequency)limit
% pubv :passband upper bound value
% plbll:passband lower bound lower (frequency)limit
% plbul:passband lower bound upper (frequency)limit
% plbv :passband lower bound value
% sll  :stopband lower (frequency)limit
% sul  :stopband upper (frequency)limit
% sv   :stopband value

maskupper = zeros(1,length(freq));
maskupper(freq<=sll | freq>=sul) = sv;
maskupper(freq>sll & freq<sul) = pubv;

masklower = zeros(1,sum(freq>plbll & freq<plbul));
masklower(:) = plbv;

end