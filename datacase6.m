function casedata = datacase6()
casedata.T = 24;
casedata.busnum = 6;
casedata.N_gen = 3;
casedata.N_ren = 2;
casedata.N_stor = 1;
% casedata.N_stor = 1;
casedata.load_pos = [4 5 6];
casedata.load_total = [4200,3960,3480,3300,3200,3600,4200,4680,4920,5280,5340,5040,4800,4560,5280,5400,5100,5340,5640,5880,6000,5400,5220,4920]/20;
load_portion = [0.3,0.3,0.4];
casedata.PDrep = load_portion'*casedata.load_total;
casedata.PDmax = 1.05 * casedata.PDrep;
casedata.PDmin = 0.95 * casedata.PDrep;
casedata.ren_pos = [3,4];
casedata.stor_pos = [5];
% casedata.stor_pos=[6,22,32,34,35];
casedata.N_load = size(casedata.load_pos,2);
casedata.N_branch = 7;
casedata.Fmax = 200*ones(casedata.N_branch,1);
casedata.Fmin = -casedata.Fmax;
PV=1*[0,0,0,0,0,0,0.0198,0.171,0.3846,0.5624,0.6842,0.7352,0.7335,0.5829,0.4298,0.2558,0.1899,0.0463,0,0,0,0,0,0]; %1MW×°»ú
WT=1*[0.709000000000000,0.684120000000000,0.630320000000000,0.581720000000000,0.569600000000000,0.551080000000000,0.476800000000000,0.437800000000000,0.462080000000000,0.465240000000000,0.467160000000000,0.459560000000000,0.456720000000000,0.436840000000000,0.427280000000000,0.445920000000000,0.410800000000000,0.362196000000000,0.398180000000000,0.481120000000000,0.506560000000000,0.546800000000000,0.640600000000000,0.700920000000000]/2.5;
casedata.ren_total = [50 50];
casedata.Rrep = zeros(casedata.N_ren,casedata.T);
casedata.Rrep(1,:) = casedata.ren_total(1,1)'*PV;
casedata.Rrep(2,:) = casedata.ren_total(1,2)'*WT;
casedata.Rmax = casedata.Rrep*1.1;
casedata.Rmin = casedata.Rrep*0.9;

casedata.stor = [
1  40  20  0.9   0.8

];
casedata.Eset = 0.5*casedata.stor(:,2)';
casedata.E0 = casedata.Eset;

casedata.startup_cost = [40,30,30];
casedata.min_up_time = [3,2,1];
casedata.min_down_time = [2,2,1];

casedata.z0 = [1,0,1];
casedata.unit_pos = [1,2,6];
casedata.pgmax = [200 180 120];
casedata.pgmin = [50 60 40];
casedata.ramp_limit = [50, 40, 40];
casedata.cost_cut=300;
casedata.cost_dcut=1000;
casedata.p0 = 0.5*casedata.pgmax .* casedata.z0;
a = [0.00056,0.00324,0.00284];
b = [8, 7.74, 8.6];
c = [150, 240, 126];
e = [100, 80, 70];
% e = [0, 0, 0];
f = [0.0419, 0.0524, 0.0785];
casedata.pcosta=zeros(3,8);
casedata.pcostb=zeros(3,8);
for i=1:casedata.N_gen
    for j = 1:8
        p1 = casedata.pgmin(i) + (j-1)*((casedata.pgmax(i)-casedata.pgmin(i))/8);
        p2 = casedata.pgmin(i) + (j)*((casedata.pgmax(i)-casedata.pgmin(i))/8);
        y1 = a(i)*p1^2 + b(i)*p1 + c(i) +e(i)*abs(sin(f(i)*(casedata.pgmin(i)-p1)));
        y2 = a(i)*p2^2 + b(i)*p2 + c(i) +e(i)*abs(sin(f(i)*(casedata.pgmin(i)-p2)));
        casedata.pcosta(i,j) = (y2-y1)/(p2-p1);
        casedata.pcostb(i,j) = y1-p1*(y2-y1)/(p2-p1);
    end
end


gama =  [0	-0.592920353982301	-0.548672566371681	-0.407079646017699	-0.477876106194690	-0.530973451327433
0	-0.407079646017699	-0.451327433628318	-0.592920353982301	-0.522123893805310	-0.469026548672566
0	0.0353982300884956	-0.743362831858407	-0.0353982300884955	-0.389380530973451	-0.654867256637168
0	0.371681415929203	0.194690265486726	-0.371681415929203	-0.0884955752212389	0.123893805309734
0	0.0353982300884954	0.256637168141593	-0.0353982300884955	-0.389380530973451	-0.654867256637168
0	-0.0353982300884956	-0.256637168141593	0.0353982300884956	-0.610619469026549	-0.345132743362832
0	-0.0353982300884956	-0.256637168141593	0.0353982300884956	0.389380530973451	-0.345132743362832];
casedata.gamad=zeros(casedata.N_branch,casedata.N_load);

for i=1:casedata.N_load
    casedata.gamad(:,i)=-gama(:,casedata.load_pos(i));
end
temp=max(abs(casedata.gamad),[],2);
temp2=find(temp<=1e-10);
gama(temp2,:)=[];
casedata.gamad(temp2,:)=[];
casedata.Fmin(temp2,:)=[];
casedata.Fmax(temp2,:)=[];
casedata.N_branch=casedata.N_branch-length(temp2);
casedata.gamag=zeros(casedata.N_branch,casedata.N_gen);
casedata.gamasto=zeros(casedata.N_branch,casedata.N_stor);
for i=1:casedata.N_gen
    casedata.gamag(:,i)=gama(:,casedata.unit_pos(i));
end
for i=1:casedata.N_stor
    casedata.gamasto(:,i)=-gama(:,casedata.stor_pos(i));
end
for i=1:casedata.N_ren
    casedata.gamar(:,i)=gama(:,casedata.ren_pos(i));
end