function [PM_u,PM_l] = phasemask_v2(type,f,val)
% This function produces a value symmetric mask dual for phase.
% type: a vector indicating the types of all subbands.
% f:    a matrix including frequency axis samples. Some rows can be smaller
% or larger than the others. Hence, smaller rows are completed with NaN
% values. These values are to be extracted before processing.
% val:  a matrix including parameters for different types. Some rows can 
% be smaller or larger than the others. Hence, smaller rows are completed 
% with NaN values. These values are to be extracted before processing.

% Produce overall frequency axis.
fo = [];
PM = [];
for i = 1:length(type)
    f_temp = f(i,:);
    f_temp(isnan(f_temp)) = [];
    val_temp = val(i,:);
    val_temp(isnan(val_temp)) = [];
    if type(i) == 0                                                        % constant
        PM = [PM,((tan(val_temp))^2)*ones(1,length(f_temp))];
    elseif type(i) == 1                                                    % linear
        PM = [PM,(tan(val_temp(1)*abs(f_temp)+val_temp(2))).^2];
    elseif type(i) == 2                                                    % quadratic %to be modified
        PM = [PM,(val_temp(1)*((f_temp-val_temp(2)).^2)).^2+val_temp(3)];
    elseif type(i) == 3                                                    % exponential %to be modified
        PM = [PM,(val_temp(1)*exp(f_temp)).^2+val_temp(2)];
    end;
    fo = [fo,f_temp];
end;
PM_u = PM;
PM_l = -PM;


% fo = [];
% PM = [];
% for i = 1:length(type)
%     f_temp = f(i,:);
%     f_temp(isnan(f_temp)) = [];
%     val_temp = val(i,:);
%     val_temp(isnan(val_temp)) = [];
%     if type(i) == 0                                                        % constant
%         PM = [PM,(((val_temp))^2)*ones(1,length(f_temp))];
%     elseif type(i) == 1                                                    % linear
%         PM = [PM,(val_temp(1)*f_temp).^2+val_temp(2)];
%     elseif type(i) == 2                                                    % quadratic %to be modified
%         PM = [PM,(val_temp(1)*((f_temp-val_temp(2)).^2)).^2+val_temp(3)];
%     elseif type(i) == 3                                                    % exponential %to be modified
%         PM = [PM,(val_temp(1)*exp(f_temp)).^2+val_temp(2)];
%     end;
%     fo = [fo,f_temp];
% end;
% PM_u = PM;
% PM_l = -PM;

end