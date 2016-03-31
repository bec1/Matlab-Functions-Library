%%
addpath('Library')
kFNlist=[];
kFFitlist=[];

%% T/T_F~0.1, good one

momimages={'03-25-2016_00_14_23_top';'03-25-2016_00_13_28_top';'03-25-2016_00_12_32_top';'03-25-2016_00_10_52_top';'03-25-2016_00_09_56_top';'03-25-2016_00_09_01_top';'03-25-2016_00_08_05_top';'03-25-2016_00_07_10_top';'03-25-2016_00_06_14_top';'03-25-2016_00_05_19_top';'03-25-2016_00_04_23_top';'03-25-2016_00_03_28_top';'03-25-2016_00_02_32_top';'03-25-2016_00_01_37_top';'03-25-2016_00_00_41_top';'03-24-2016_23_59_46_top';'03-24-2016_23_58_50_top';'03-24-2016_23_57_55_top';'03-24-2016_23_56_59_top';'03-24-2016_23_56_04_top';'03-24-2016_23_55_08_top';'03-24-2016_23_54_13_top';'03-24-2016_23_53_17_top';'03-24-2016_23_52_22_top';'03-24-2016_23_51_26_top';'03-24-2016_23_50_31_top';'03-24-2016_23_49_10_top';'03-24-2016_23_48_14_top';'03-24-2016_23_47_19_top';'03-24-2016_23_46_24_top';'03-24-2016_23_45_28_top';'03-24-2016_23_44_33_top'};
bgimages={'03-24-2016_22_24_18_top';'03-24-2016_22_23_27_top';'03-24-2016_22_22_13_top'};

output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4,'TailRange',[50,450]);
%%
kFNlist=[kFNlist,output.kF_num];
kFFitlist=[kFFitlist,output.kF_Fit];

%% T/T_F~0.5, gpood one
momimages={'03-24-2016_23_31_27_top';'03-24-2016_23_30_32_top';'03-24-2016_23_29_36_top';'03-24-2016_23_27_50_top';'03-24-2016_23_25_58_top';'03-24-2016_23_25_03_top';'03-24-2016_23_24_07_top';'03-24-2016_23_23_12_top';'03-24-2016_23_22_16_top';'03-24-2016_23_21_20_top';'03-24-2016_23_20_25_top';'03-24-2016_23_19_29_top';'03-24-2016_23_18_34_top';'03-24-2016_23_17_38_top';'03-24-2016_23_16_42_top';'03-24-2016_23_15_47_top';'03-24-2016_23_14_51_top';'03-24-2016_23_12_17_top';'03-24-2016_23_11_21_top';'03-24-2016_23_10_25_top';'03-24-2016_23_09_30_top';'03-24-2016_23_08_34_top';'03-24-2016_23_07_39_top';'03-24-2016_23_06_43_top';'03-24-2016_23_05_47_top';'03-24-2016_23_04_52_top';'03-24-2016_23_02_35_top';'03-24-2016_23_01_40_top';'03-24-2016_23_00_44_top';'03-24-2016_22_59_48_top'};

bgimages={'03-24-2016_23_26_54_top';'03-24-2016_22_24_18_top';'03-24-2016_22_23_27_top';'03-24-2016_22_22_13_top'};

output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4,'TailRange',[40,420]);

%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]

%% T/T_F~0.17
momimages={'03-24-2016_21_06_46_top';'03-24-2016_21_05_55_top';'03-24-2016_21_05_03_top';'03-24-2016_21_04_12_top';'03-24-2016_21_03_20_top';'03-24-2016_21_02_29_top';'03-24-2016_21_01_38_top';'03-24-2016_21_00_00_top';'03-24-2016_20_59_09_top';'03-24-2016_20_58_18_top';'03-24-2016_20_55_44_top';'03-24-2016_20_54_52_top';'03-24-2016_20_53_09_top';'03-24-2016_20_52_18_top';'03-24-2016_20_51_27_top';'03-24-2016_20_50_35_top';'03-24-2016_20_49_44_top';'03-24-2016_20_48_52_top';'03-24-2016_20_48_01_top'};
bgimages={'03-24-2016_21_00_00_top';'03-24-2016_20_59_09_top';'03-24-2016_20_58_18_top'};
output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4,'TailRange',[40,420]);
%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]

%% Fitting not working well T/T_F~0.2
momimages={'03-24-2016_22_40_59_top';'03-24-2016_22_37_32_top';'03-24-2016_22_36_41_top';'03-24-2016_22_35_50_top';'03-24-2016_22_34_59_top';'03-24-2016_22_34_08_top';'03-24-2016_22_32_39_top';'03-24-2016_22_31_48_top';'03-24-2016_22_30_57_top';'03-24-2016_22_30_06_top';'03-24-2016_22_28_44_top';'03-24-2016_22_27_53_top';'03-24-2016_22_27_02_top';'03-24-2016_22_24_18_top';'03-24-2016_22_23_27_top';'03-24-2016_22_22_13_top';'03-24-2016_22_21_22_top';'03-24-2016_22_20_31_top';'03-24-2016_22_18_34_top';'03-24-2016_22_17_43_top';'03-24-2016_22_16_51_top';'03-24-2016_22_15_46_top';'03-24-2016_22_14_43_top';'03-24-2016_22_13_23_top';'03-24-2016_22_12_19_top'};
bgimages={'03-24-2016_21_00_00_top';'03-24-2016_20_59_09_top';'03-24-2016_20_58_18_top'};
output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4);
%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]

%% T/T_F~0.2, semi-goodone

momimages={'03-24-2016_21_29_23_top';'03-24-2016_21_28_32_top';'03-24-2016_21_27_40_top';'03-24-2016_21_26_49_top';'03-24-2016_21_25_57_top';'03-24-2016_21_25_06_top';'03-24-2016_21_24_14_top';'03-24-2016_21_23_23_top';'03-24-2016_21_22_31_top';'03-24-2016_21_21_40_top';'03-24-2016_21_20_48_top';'03-24-2016_21_19_57_top'};
bgimages={'03-24-2016_21_00_00_top';'03-24-2016_20_59_09_top';'03-24-2016_20_58_18_top'};
output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4);

%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]

%% T/T_F~0.17
momimages={'03-22-2016_00_44_20_top';'03-22-2016_00_42_07_top';'03-22-2016_00_41_12_top';'03-22-2016_00_40_18_top';'03-22-2016_00_39_24_top';'03-22-2016_00_38_30_top';'03-22-2016_00_37_35_top';'03-22-2016_00_36_41_top';'03-22-2016_00_35_46_top';'03-22-2016_00_34_15_top'};
bgimages={'03-21-2016_23_40_31_top'}; 
output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4);
%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]

%% T/T_F~0.15
momimages={'03-19-2016_03_03_48_top';'03-19-2016_03_02_53_top';'03-19-2016_03_01_59_top';'03-19-2016_03_01_05_top';'03-19-2016_03_00_10_top';'03-19-2016_02_59_16_top';'03-19-2016_02_58_22_top';'03-19-2016_02_57_28_top';'03-19-2016_02_56_33_top';'03-19-2016_02_55_39_top';'03-19-2016_02_54_44_top';'03-19-2016_02_53_50_top';'03-19-2016_02_52_56_top';'03-19-2016_02_52_02_top';'03-19-2016_02_51_07_top';'03-19-2016_02_50_13_top';'03-19-2016_02_49_19_top';'03-19-2016_02_48_25_top';'03-19-2016_02_47_30_top'};
bgimages={'03-19-2016_02_39_33_top';'03-19-2016_02_38_39_top';'03-19-2016_02_37_45_top';'03-19-2016_02_36_39_top'};
output=momentumfocusRV(momimages,bgimages,'Nbins' ,60,'SM',4);
%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]

%% T/T_F~0.1, goodone
momimages={'03-12-2016_01_54_06_top';'03-12-2016_01_53_18_top';'03-12-2016_01_52_31_top';'03-12-2016_01_51_43_top';'03-12-2016_01_49_03_top';'03-12-2016_01_48_16_top';'03-12-2016_01_46_41_top';'03-12-2016_01_45_07_top';'03-12-2016_01_40_40_top'};
bgimages={'03-12-2016_01_19_51_top'};
output=momentumfocusRV(momimages,bgimages,'Nbins' ,30,'SM',2);

%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]
%% T/T_F~0.17
momimages={'03-04-2016_21_31_32_top';'03-04-2016_21_30_39_top';'03-04-2016_21_29_45_top';'03-04-2016_21_28_51_top';'03-04-2016_21_27_57_top';'03-04-2016_21_27_03_top';'03-04-2016_21_26_10_top';'03-04-2016_21_25_16_top';'03-04-2016_21_24_22_top';'03-04-2016_21_18_09_top';'03-04-2016_21_17_16_top';'03-04-2016_21_16_22_top';'03-04-2016_21_15_28_top';'03-04-2016_21_14_34_top';'03-04-2016_21_13_40_top';'03-04-2016_21_12_47_top';'03-04-2016_21_11_53_top';'03-04-2016_21_10_59_top';'03-04-2016_21_10_05_top';'03-04-2016_21_09_12_top';'03-04-2016_21_08_18_top';'03-04-2016_21_07_24_top';'03-04-2016_21_06_30_top';'03-04-2016_21_05_37_top';'03-04-2016_21_04_43_top';'03-04-2016_21_02_55_top';'03-04-2016_21_02_01_top'};
bgimages={'03-04-2016_21_22_44_top';'03-04-2016_21_21_28_top';'03-04-2016_21_20_34_top';'03-04-2016_21_19_41_top';'03-04-2016_21_36_16_top';'03-04-2016_21_35_23_top';'03-04-2016_21_34_29_top';'03-04-2016_21_33_35_top';'03-04-2016_21_43_27_top';'03-04-2016_21_42_33_top';'03-04-2016_21_41_39_top';'03-04-2016_21_40_45_top';'03-04-2016_21_39_51_top';'03-04-2016_21_38_58_top'};
output=momentumfocusRV(momimages,bgimages,'Nbins' ,50,'SM',3);

%%
kFNlist=[kFNlist,output.kF_num]
kFFitlist=[kFFitlist,output.kF_Fit]


%%
scatter(kFNlist,kFFitlist,'DisplayName','Experimental Data');
hold on
line([0,4e6],[0,4e6],'DisplayName','y=x','color','r','linewidth',2);
xlabel('k_F from total atom number (m^{-1})');
ylabel('k_F from Fermi-Dirac distribution (m^{-1})');
legend('show')
hold off