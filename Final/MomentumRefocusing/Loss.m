addpath('Library');
imagenames={'03-11-2016_20_34_36_top';'03-11-2016_20_33_49_top';'03-11-2016_20_33_01_top';'03-11-2016_20_32_14_top';'03-11-2016_20_31_26_top';'03-11-2016_20_30_39_top';'03-11-2016_20_27_46_top';'03-11-2016_20_26_59_top';'03-11-2016_20_26_12_top';'03-11-2016_20_25_25_top';'03-11-2016_20_24_37_top';'03-11-2016_20_23_50_top';'03-11-2016_20_23_03_top';'03-11-2016_20_19_39_top';'03-11-2016_20_18_52_top';'03-11-2016_20_18_05_top';'03-11-2016_20_17_17_top';'03-11-2016_20_16_30_top';'03-11-2016_20_15_43_top';'03-11-2016_20_14_56_top';'03-11-2016_20_14_08_top';'03-11-2016_20_13_21_top';'03-11-2016_20_12_33_top';'03-11-2016_20_11_46_top';'03-11-2016_20_10_59_top';'03-11-2016_20_10_00_top';'03-11-2016_20_08_42_top';'03-11-2016_20_07_55_top';'03-11-2016_20_05_52_top';'03-11-2016_20_05_05_top';'03-11-2016_20_04_17_top'};
time=[30;29;28;27;26;25;24;23;22;21;20;19;18;17;16;15;14;13;12;11;10;9;8;7;6;5;4;3;2;1;0];
% imagenames={'03-12-2016_01_33_02_top';'03-12-2016_01_30_22_top';'03-12-2016_01_29_34_top';'03-12-2016_01_28_47_top';'03-12-2016_01_27_59_top';'03-12-2016_01_27_12_top';'03-12-2016_01_26_24_top';'03-12-2016_01_25_37_top';'03-12-2016_01_24_49_top';'03-12-2016_01_24_02_top';'03-12-2016_01_23_15_top';'03-12-2016_01_22_27_top';'03-12-2016_01_21_40_top'};
% time=[26;24;22;20;18;16;14;12;10;8;4;2;0.00100000000000000];
Nsat=330;
pixellength=1.44*10^(-6);
sigma0=0.215*10^(-12)/2;
ROI=[200,2,130,500];
%%
N=length(imagenames);
num=time*0;

for i=1:N
    [~,tempraw]=imagedata(imagenames{i});
    Ntemp=AtomNumber(tempraw,pixellength.^2,sigma0, Nsat);
    Ntemp=imcrop(Ntemp,ROI);
    n=sum(Ntemp,2)';
    z=1:length(n);
%     h=figure();
%     scatter(z,n);
%     questdlg('Now give the range for tail fitting');
%     [x,y]=getpts(h);
%     close(h);
    zmin=50;
    zmax=390;
    n=TailTailor(n,z,zmin,zmax);
    num(i)=sum(n);
end

%%
scatter(time,num)