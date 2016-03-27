addpath('Library');
bgimages1={'03-22-2016_00_50_45_top';'03-22-2016_00_43_01_top';'03-21-2016_23_40_31_top'};
momimages1={'03-22-2016_01_00_49_top';'03-22-2016_00_59_54_top';'03-22-2016_00_58_23_top';'03-22-2016_00_57_28_top';'03-22-2016_00_56_10_top';'03-22-2016_00_55_16_top';'03-22-2016_00_54_22_top';'03-22-2016_00_53_27_top';'03-22-2016_00_52_33_top';'03-22-2016_00_51_39_top';'03-22-2016_00_49_50_top';'03-22-2016_00_48_56_top'};
%bgimages1={'03-19-2016_02_39_33_top';'03-19-2016_02_38_39_top';'03-19-2016_02_37_45_top';'03-19-2016_02_36_39_top'};
%momimages1={'03-19-2016_03_03_48_top';'03-19-2016_03_02_53_top';'03-19-2016_03_01_59_top';'03-19-2016_03_01_05_top';'03-19-2016_03_00_10_top';'03-19-2016_02_59_16_top';'03-19-2016_02_58_22_top';'03-19-2016_02_57_28_top';'03-19-2016_02_56_33_top';'03-19-2016_02_55_39_top';'03-19-2016_02_54_44_top';'03-19-2016_02_53_50_top';'03-19-2016_02_52_56_top';'03-19-2016_02_52_02_top';'03-19-2016_02_51_07_top';'03-19-2016_02_50_13_top';'03-19-2016_02_49_19_top';'03-19-2016_02_48_25_top';'03-19-2016_02_47_30_top'};
%%
output=momentumfocusRV(momimages1,bgimages1,'Nbins' ,50,'SM',2);
n=output{1};
kz=output{2}/1e6;

