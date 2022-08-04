function [csi] = csi_gen(Nu,Nsc,var,ds,d0)
%Channel State Information Generateror Fuction
%Generates random CSI for a given number of users and subcarriers
%Time varying fading is according to a first order auto-regressive moodel

%DATA in needed to perform the technique
% ----------------------------------------------------------------------
%--Nu is the number of users in the system
%--Nsc is the number of subcarriers in the system
%--var is the csi variance.
%--ds is the distance moved between each csi sample
%--d0 is the decay factor
% ----------------------------------------------------------------------


%DATA out
% ----------------------------------------------------------------------
%--csi cell data of csi for all users
% ----------------------------------------------------------------------
    csi = cell(Nsc,Nu);
    alpha = exp(-(ds/d0));
    fade_var=(1-alpha^2)*var;
    for i=1:Nu
        temp = ((abs(normrnd(0,var,[1,Nsc])).^2)/var);
        csi{1,i}(1:Nsc) = temp;
        for j=2:Nsc
            %First order autoregressive time varying fading channel
            csi{j,i}(1:Nsc) = alpha*csi{j-1,i}(1:Nsc)+normrnd(0,fade_var,[1,Nsc]);
            %Don't allow csi to be less than 0.
            csi{j,i}(1:Nsc) = max(csi{j,i}(1:Nsc),0);
        end
    end
end

