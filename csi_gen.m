function [csi] = csi_gen(Nu,Nsc,var,alpha,fade_var)
%Channel State Information Generateror Fuction
%Generates random CSI for a given number of users and subcarriers

%DATA in needed to perform the technique
% ----------------------------------------------------------------------
%--Nu is the number of users in the system
%--Nsc is the number of subcarriers in the system
%--var is the csi variance.
% ----------------------------------------------------------------------


%DATA out
% ----------------------------------------------------------------------
%--csi cell data of csi for all users
% ----------------------------------------------------------------------
    csi = cell(Nsc,Nu);
    for i=1:Nu
        temp = ((abs(normrnd(0,var,[1,Nsc])).^2)/var);
        csi{1,i}(1:Nsc) = temp;
        for j=2:Nsc
            %First order autoregressive time varying fading channel
            csi{j,i}(1:Nsc) = alpha*csi{j-1,i}(1:Nsc)+normrnd(0,fade_var,[1,Nsc]);
        end
    end
end

