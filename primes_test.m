clear;
runs = primes(15);
%runs = [0,17];
output = [];
d0 = 0.1;
ds = 0.5;
var = 1.5;
for r=runs(2:end)
    s.Rb_size=12;
    s.Nu = r;
    p_aloc_rate = s.Nu;
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
    [s.L,s.C,s.L_nh,s.C_nh,s.C_ideal] = chaoticmap(s.Nu,s.Nsc,s.RBs,s.alpha,s.rx,s.ry,s.Kc,s.Rb_size,s.csi);
    [s.Pn_opt,s.csi_ra,s.Cn] = waterfilling(s.csi,s.C,s.RBs,s.Pt,s.Rb_size,p_aloc_rate);
    [s.Pn_opt2,s.csi_ra2,s.Cn2] = waterfilling(s.csi,s.C_nh,s.RBs,s.Pt,s.Rb_size,p_aloc_rate);
    [s.Pn_opt3,s.csi_ra3,s.Cn3] = waterfilling(s.csi,s.C_ideal,s.RBs,s.Pt,s.Rb_size,p_aloc_rate);
    [s.Pn_opt4,s.csi_ra4,s.Cn4] = waterfilling(s.csi,s.C_ideal,s.RBs,s.Pt,s.Rb_size,1);

    s.count = zeros(1,s.Nu);
    for i=1:s.Nu
        s.count(i)=sum(sum(s.C==i));
        %fprintf("user %d clusters %d\n", i, s.count(i));
    end

    %use the throughput for user 1
    %s.Th=sum(cell2mat(s.Cn))/s.Nu/s.Nclu;
    s.Th=sum(s.Cn{1})/s.Nclu;
    fprintf("FH average throughput per slot %f bps/Hz\n", s.Th);
    %s.Th2=sum(cell2mat(s.Cn2))/s.Nu/s.Nclu;
    s.Th2=sum(s.Cn2{1})/s.Nclu;
    s.Nts=size(s.C,2);
    fprintf("Non FH average throughput per slot %f bps/Hz\n", s.Th2);
    %s.Th3=sum(cell2mat(s.Cn3))/s.Nu/s.Nclu;
    s.Th3=sum(s.Cn3{1})/s.Nclu;
    fprintf("Ideal FH average throughput per slot %f bps/Hz\n", s.Th3)
    s.Th4=sum(s.Cn4{1})/s.Nclu;
    fprintf("Ideal FH Continous Power average throughput per slot %f bps/Hz\n", s.Th4)
    output = [output s];
end
%save('newstruct.mat', 'output')

%Plotting results for first time slot
f1 = figure(1);
clf(f1);
%hold on;
sz=size(s.C);
%sz=[4,5]
%for i=1:sz(2)
%    scatter(C(:,i),1:1:size(C,1),'filled','DisplayName',['timeslot ',num2str(i)])
%end
heatmap(s.C);
title('Resource Allocation for each User Freqency Hopping')
ylabel('Resource cluster')
xlabel('Timeslot')
%legend();

f2 = figure(2);
clf(f2);
bar(s.Pn_opt{1,1},1,'r')
xlabel('Subchannels')
ylabel('Power Allocated')
title('Power Allocation. Total Power for the channel = ',num2str(sum(s.Pn_opt{1,1})))

f3 = figure(3);
clf(f3);
bar(s.csi_ra{1,1},1)
xlabel('Subchannels')
ylabel('CSI')
title('CSI-RA for the first Time Slot')

f4 = figure(4);
clf(f4);
%hold on;
sz=size(s.C_nh);
%sz=[4,5]
%for i=1:sz(2)
%    scatter(C(:,i),1:1:size(C,1),'filled','DisplayName',['timeslot ',num2str(i)])
%end
heatmap(s.C_nh);
title('Resource Allocation for each User No Hopping')
ylabel('Resource cluster')
xlabel('Timeslot')
%legend();

f5 = figure(5);
clf(f5);
%hold on;
sz=size(s.C_ideal);
heatmap(s.C_ideal);
title('Resource Allocation for each User Ideal')
ylabel('Resource cluster')
xlabel('Timeslot')

%f1 = figure(1);
%clf(f1);
% set(f1, 'Color', [1 1 0.1])
%bar(s.C(:,1),1,'g')
%title('Resource Allocation for the first Symbol Period')
%xlabel('Subchannels/Subcarriers')
%ylabel('Users')