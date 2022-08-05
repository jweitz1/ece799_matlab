function loss_test(d0, N)
    d0 = str2num(d0);
    N = str2num(N);
    output = [];
    output2 = [];
    output3 = [];
    output4 = [];
    th1 = [];
    th2 = [];
    th3 = [];
    th4 = [];
    th5 = [];
    th6 = [];
    ds_array = linspace(0,1,25);
    var = 1.5;
    k = 50; % inner loop count
    for ds = ds_array
        percentage = 0;
        percentage2 = 0;
        percentage3 = 0;
        percentage4 = 0;
        t1 = 0;
        t2 = 0;
        t3 = 0;
        t4 = 0;
        t5 = 0;
        t6 = 0;
        for i = 1:k
            s.Rb_size=12;
            s.Nu = N;
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
            [s.Pn_opt5,s.csi_ra5,s.Cn5] = waterfilling(s.csi,s.C_nh,s.RBs,s.Pt,s.Rb_size,1);
            [s.Pn_opt6,s.csi_ra6,s.Cn6] = waterfilling(s.csi,s.C,s.RBs,s.Pt,s.Rb_size,1);

            s.count = zeros(1,s.Nu);
            for i=1:s.Nu
                s.count(i)=sum(sum(s.C==i));
                %fprintf("user %d clusters %d\n", i, s.count(i));
            end

            %use the throughput for user 1
            s.Th=sum(cell2mat(s.Cn))/s.Nu/s.Nclu;
            fprintf("FH average throughput per slot %f bps/Hz\n", s.Th);
            s.Th2=sum(cell2mat(s.Cn2))/s.Nu/s.Nclu;
            s.Nts=size(s.C,2);
            fprintf("Non FH average throughput per slot %f bps/Hz\n", s.Th2);
            s.Th3=sum(cell2mat(s.Cn3))/s.Nu/s.Nclu;
            fprintf("Ideal FH average throughput per slot %f bps/Hz\n", s.Th3)
            s.Th4=sum(cell2mat(s.Cn4))/s.Nu/s.Nclu;
            fprintf("Ideal FH Continous Power average throughput per slot %f bps/Hz\n", s.Th4)
            s.Th5=sum(cell2mat(s.Cn5))/s.Nu/s.Nclu;
            fprintf("Ideal FH Continous Power average throughput per slot %f bps/Hz\n", s.Th5)
            s.Th6=sum(cell2mat(s.Cn6))/s.Nu/s.Nclu;
            fprintf("FH Continous Power average throughput per slot %f bps/Hz\n", s.Th6)
            %fprintf("Hopping Throughput Loss %f%%\n", 100.0-(s.Th/s.Th2)*100.0);
            percentage = percentage + (100.0-(abs(s.Th)/abs(s.Th2))*100.0);
            percentage2 = percentage2 + (100.0-(abs(s.Th)/abs(s.Th3))*100.0);
            percentage3 = percentage3 + (100.0-(abs(s.Th)/abs(s.Th4))*100.0);
            percentage4 = percentage4 + (100.0-(abs(s.Th)/abs(s.Th5))*100.0);
            %output = [output s];
            t1 = t1 + s.Th;
            t2 = t2 + s.Th2;
            t3 = t3 + s.Th3;
            t4 = t4 + s.Th4;
            t5 = t5 + s.Th5;
            t6 = t6 + s.Th6;
        end
        percentage = percentage/k;
        percentage2 = percentage2/k;
        percentage3 = percentage3/k;
        percentage4 = percentage4/k;
        t1 = t1/k;
        t2 = t2/k;
        t3 = t3/k;
        t4 = t4/k;
        t5 = t5/k;
        t6 = t6/k;
        fprintf("After %d runs, for ds = %f:\n", k, ds);
        fprintf("CSM vs Non Hopping Throughput Loss %f%%\n", percentage);
        fprintf("CSM vs Ideal Hopping Throughput Loss %f%%\n", percentage2);
        fprintf("CSM vs Ideal Hopping Continuous Power Allocation Throughput Loss  %f%%\n", percentage3);
        fprintf("CSM vs Non Hopping Continuous Power Allocation Throughput Loss  %f%%\n", percentage4);
        output = [output percentage];
        output2 = [output2 percentage2];
        output3 = [output3 percentage3];
        output4 = [output4 percentage4];
        th1 = [th1 t1];
        th2 = [th2 t2];
        th3 = [th3 t3];
        th4 = [th4 t4];
        th5 = [th5 t5];
        th6 = [th6 t6];
        f = sprintf('loss_comp_N%d_d0_0p%d.mat', N, round(10*abs(d0 - fix(d0))));
        save(f, 'd0', 'ds_array', 'output', 'output2', 'output3', 'output4','N', 'k', 'th1', 'th2', 'th3', 'th4', 'th5', 'th6')
    end
    %figure(1);
    %plot(ds_array, output);
    %title('Hopping Throughput Loss vs d_{s}');
    %ylabel('Througput Loss %');
    %xlabel('d_{s}');
end
