function [L,C,csi,L_nh,C_nh] = chaoticmap(Nu,Nsc,RBs,alpha,rx,ry,Kc,Rb_size, csi)

 %When Frequency Hopping is activated, this code is executed 
% after run button is pushed. Currently, Chaotic Map Technique 
% is used but other techniques can be implemented in future
% It would be interesting to implement a list of techniques 
% to make the user able to select the one desired.
%We should create a listBox with a set of FH techniques
% Then execute the FH code corresponding with the FH technique selected. 


%Frequency Hopping Technique based on chaotic map technique

%DATA in needed to perform the technique
% ----------------------------------------------------------------------
%--kusers is the number of users in the system
%--Nsc is the number of subcarriers in the system
%--Rbs is the number of Resource Blocks for each user
%--alpha represents a pattern index for alpha=2,...,Nclu-1. alpha=2 is an
%  example used for simulations
%--rx, ry are integers varying from 0 to Nclu-1.(rx,ry)=(3,1) are examples
%  used for simulations
%--Kc is a positive integer. 1000 is an example used for simulations.
% ----------------------------------------------------------------------


%DATA out
% ----------------------------------------------------------------------
%--CSM_all contains the FH pattern for each user K, it means the Frequency
%       hops from the subcarriers in matrix L of user K to the subcarriers in
%       matrix C
%--L represents the first FH assignation to the users, each value of the
%       matrix represents the assignation of Nclu clusters
%       to the user marked for the matrix. Seed matrix.
%--C represents the change in the FH assignation based on chaotic map
% ----------------------------------------------------------------------

    Nclu=floor(Nsc/(RBs*Rb_size)); %Number of clusters in the channel
                              %It is RBs * Rb_size because each RB has Rb_size
                              %subcarriers
                              
    %csi = cell(1,Nu);
    %for i=1:Nu
    %    %csi{1,i}(1:Nsc) = rand(1,Nsc);
    %    %distance =  normrnd(1.0,0.5);
    %    %temp = ((abs(normrnd(1,0.5,[1,Nsc])).^2)./0.5)*distance;
    %    var = 1.5;
    %    temp = ((abs(normrnd(0,var,[1,Nsc])).^2)/var);
    %    %temp(temp < 0)=0.000001;
    %    %temp(temp > 1)=1;
    %    csi{1,i}(1:Nsc) = temp;
    %end

    %Next code will be executed only if FH is activated

    %if useHop==1

        L=zeros(Nclu,Nclu);             %Seed Matrix

        %Compute alpha inverse
        alpha_inv=0;    
        for n=0:Nclu-1
            if mod(alpha*n,Nclu)==1
                alpha_inv=n;
            end
        end

        %Construct an Nclu x Nclu seed matrix using Latin
        %Square
        for i = 1:Nclu
            for j=1:Nclu
                L(i,j)=mod((alpha*(i-1)+(j-1)),Nclu);
            end
        end
        kusers=max(max(L))+1;                    %Max users in the system
        S_FH_all=cell(1,kusers);                 %FH pattern matrix for each user
        s_fh_all=cell(1,kusers);                 %FH pattern matrix for each user
        CSM_all=cell(1,kusers);                  %CSM FH pattern matrix for each user

        %S_FH represents the positions in L for user k
        for k=1:kusers                 
            row=1;
            column=1;
            for i = 1:Nclu
                for j=1:Nclu
                    if L(i,j)+1==(k)
                        S_FH(row,column)=i;
                        S_FH(row,column+1)=j;
                        column=1;
                        row=row+1;
                    end
                end
            end

            s_fh=zeros(length(S_FH),1);
            for j=0:Nclu-1
                s_fh(j+1)=(mod(floor(((k-1)-(j))*alpha_inv), Nclu));
            end
            s_fh=s_fh';

            %Convert the FH pattern into CSM FH pattern
            %Rows CSM Matrix
            for j=0:Nclu-1
                CSM_row(j+1)=(mod(s_fh(j+1)+j+rx+ry,Nclu));
            end                   
            %Columns CSM Matrix
            for j=0:Nclu-1
                CSM_column(j+1)=(mod(j+ry+round((Kc*sin((CSM_row(j+1)*Nclu)/(2*pi)))),Nclu));    
            end 
            %CSM Matrix
            CSM=zeros(size(S_FH,1),size(S_FH,2));
            CSM(:,1)=CSM_row(:)+1; %Scaled to +1 because matlab does not have 0 indexes
            CSM(:,2)=CSM_column(:)+1;

            S_FH_all(1,k)={S_FH};
            s_fh_all(1,k)={s_fh};
            CSM_all(1,k)={CSM};
        end

        %Seed Matrix converted according to the FH pattern
        C=L;
        for i=1:length(CSM_all)
            for j=1:size(S_FH_all{1,i},1) 
                C(CSM_all{1,i}(j,1),CSM_all{1,i}(j,2))=L(S_FH_all{1,i}(j,1),S_FH_all{1,i}(j,2));     
            end
        end

        %If Nusers is lower than Nclusters, we randomly assign 
        %the free subcarriers between the users. 
        if Nu<Nclu
            for i=1:size(C,1)
                for j=1:size(C,2)
                    if Nu<=C(i,j)
                        C(i,j)=randi([1 Nu],1);                               
                    end
                end
            end
        %We cannot have more users than clusters
        elseif Nu>Nclu
           error('Max Users allowed = #Subcarriers/(#RB*Rb_size)');
        end

        %Not zeros within the matrix (user one instead of user zero) 
        C=C+1;


    %If not FH is desired
    %else

        C_nh=zeros(Nclu,Nclu);
        C_nh=C_nh+3301;
        %First RA approach when no FH is activated
        %It would be interesting to compare the CSI between all
        %the users but it requires a high computation cost when
        %there are a lot of subcarriers.
        cluster_avail = ones(1,Nclu);
        for column=1:Nu
            cluster_qual = zeros(1,Nclu);
            for row=1:Nclu
                cluster_qual(row) = sum(csi{1,column}((row-1)*RBs*Rb_size+1:row*RBs*Rb_size));
                %if sum(csi{1,1}((row-1)*RBs*Rb_size+1:row*RBs*Rb_size))>sum(csi{1,column}((row-1)*RBs*Rb_size+1:row*RBs*Rb_size))
                %    C_nh(row,column)=0;
                %else
                %    C_nh(row,column)=column-1;
                %end
            end
            cluster_qual(cluster_avail==0)=0;
            [M, I] = max(cluster_qual);
            C_nh(I,:)=column;
            cluster_avail(I) = 0;
            
        end


        if Nu<Nclu
            for i=1:size(C_nh,1)
                for j=1:size(C_nh,2)
                    if Nu<=C_nh(i,j)
                        C_nh(i,j)=randi([1 Nu],1);                               
                    end
                end
            end

        elseif Nu>Nclu
           error('Max Users allowed = #Subcarriers/(#RB*Rb_size)');
        end

        %C_nh=C_nh+1;
        L_nh=C_nh;
    %end
    %When there are less users than clusters we assign randomly
    %the free clusters to the users
    if Nu<Nclu
        for i=1:size(C,1)
            for j=1:size(C,2)
                if Nu<=C(i,j)
                    C(i,j)=randi([1 Nu],1);                               
                end
            end
        end

    elseif Nu>Nclu
        error('Max Users = #Subcarriers/(#RB*Rb_size)');
    end
    

end

