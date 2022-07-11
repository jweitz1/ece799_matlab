clear;
N = 17;
output = [];
d0 = 0.5;
ds_array = linspace(0,1,25);
var = 1.5;
k = 50; % inner loop count
for ds = ds_array
    percentage = 0;
    for i = 1:k
        s.Rb_size=12;
        s.Nu = N;
        s.Nsc = s.Nu*s.Rb_size;
        s.RBs = floor(s.Nsc/s.Rb_size/s.Nu);
        s.Nclu=floor(s.Nsc/(s.RBs*s.Rb_size));
        s.alpha = randi([2,s.Nclu-1],1);
        s.rx = randi([0,s.Nclu-1],1);
        s.ry = randi([0,s.Nclu-1],1);
        s.Kc = randi([0,10000000],1);
        s.Pt=1000;
        fprintf("RBs %d Nclu %d alpha %d rx %d ry %d Kc %d\n", s.RBs, s.Nclu, s.alpha, s.rx, s.ry, s.Kc);
        s.csi = csi_gen(s.Nu,s.Nsc,var,ds,d0);
        [s.L,s.C,s.L_nh,s.C_nh] = chaoticmap(s.Nu,s.Nsc,s.RBs,s.alpha,s.rx,s.ry,s.Kc,s.Rb_size,s.csi);
        [s.Pn_opt,s.csi_ra,s.Cn] = waterfilling(s.csi,s.C,s.RBs,s.Pt,s.Rb_size);
        [s.Pn_opt2,s.csi_ra2,s.Cn2] = waterfilling(s.csi,s.C_nh,s.RBs,s.Pt,s.Rb_size);

        s.count = zeros(1,s.Nu);
        for i=1:s.Nu
            s.count(i)=sum(sum(s.C==i));
            %fprintf("user %d clusters %d\n", i, s.count(i));
        end

        s.Th=sum(cell2mat(s.Cn))/s.Nu/s.Nclu;
        fprintf("FH Average throughput per slot %f bps/Hz\n", s.Th);
        s.Th2=sum(cell2mat(s.Cn2))/s.Nu/s.Nclu;
        s.Nts=size(s.C,2);
        fprintf("Non FH Average throughput per slot %f bps/Hz\n", s.Th2);
        %fprintf("Hopping Throughput Loss %f%%\n", 100.0-(s.Th/s.Th2)*100.0);
        percentage = percentage + (100.0-(abs(s.Th)/abs(s.Th2))*100.0);
        %output = [output s];
    end
    percentage = percentage/k;
    fprintf("After %d runs, for ds = %f:\n", k, ds);
    fprintf("Hopping Throughput Loss %f%%\n", percentage);
    output = [output percentage];
end
figure(1);
plot(ds_array, output);
title('Hopping Throughput Loss vs d_{s}');
ylabel('Througput Loss %');
xlabel('d_{s}');
%save('newstruct.mat', 'output')