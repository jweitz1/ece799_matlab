function [Pn_opt,csi_ra, Cn] = waterfilling(csi,C, RBs, Pt, Rb_size, update_period)

    % Water Filling Algorithm
    % ---------------------------------------------------------------

    % - Water Filling algorithm to find the optimal value to allocate power for
    % different channels. 
    % - CSI (Subchannel information) is assumed to be known. csi is
    % a cell containing CSI information for all the users in each
    % subcarrier
    % - C is the Matrix with the allocation for each user
    % - RBs is the number of resoruce blocks per user
    % - Pt is total power in the system
    % - Maximum Throughput, Cn, for each Time Slot
    % - In future, Bandwidth can be included (B), where B represents the bandwidth used for
    % subchannels.

    % --------------------------------------------------------------

    % Instead of using loops, we are going to use matrix vector computation because
    % MATLAB is very slow using loops but very fast when using matrix vector
    % computation.


    %CSI of each user is computed for the RA selected.            
    Nclu = size(C,1);

    csi_ra=cell(Nclu,Nclu); %CSI for the current RA
    for column=1:Nclu
        for row=1:Nclu
            for user=1:Nclu
                csi_ra{user,column}((row-1)*RBs*Rb_size+1:row*RBs*Rb_size) = csi{user,C(row,column)}((row-1)*RBs*Rb_size+1:row*RBs*Rb_size);                
            end
        end
    end

    %mc = length(csi); % Number of subchannels/subcarriers
    M = 2e3;          % Number of grid points you want to compute  the lagrangian dual g of mu
    mu_axis= linspace(1e-15,max(max(cell2mat(csi))),M);

    % Power allocation for each Time Slot
    Pn = cell(1,Nclu);          % Power
    g = cell(1,Nclu);           % Cost Function
    Pn_opt = cell(1,Nclu);      % Optimal power
    Cn = cell(1,Nclu);          % Maximum Rate reached

    %For each user
    for i=1:Nclu
        Cn{1,i}=0;
        %For each time slot
        for slot=1:Nclu         
            if mod(slot-1,update_period)==0
                %If it's the update period, update the power allocation for the slot
                Pn{1,i} = max((1./mu_axis - (1./csi_ra{i,slot})'),0); % A value less than 0 will be replaced by 0.
                %% Lagrange Dual Function
                g{1,i} = sum(log2(1+Pn{1,i}.*(repmat(csi_ra{i,slot}',[1 M])))) - mu_axis.*(sum(Pn{1,i}) - Pt);
                %We have to find the minimum of g
                [ind]= find(g{1,i}==min(g{1,i})); %We need the index to compute the optimal mu
                mu = mu_axis(ind);
                %Compute the optimal Powers based on the optimal mu value.
                Pn_opt{1,i} = max(1./mu - 1./csi_ra{i,slot},0);
            end
            
            Cn{1,i} = Cn{1,i}+sum(1*log2(1+Pn_opt{1,i}.*csi_ra{i,slot})); % Max Throughput without considering the bandwidth, B
        end
    end
end
