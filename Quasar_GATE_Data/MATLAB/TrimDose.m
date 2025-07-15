%the MHD file was imported to Python, where it was converted to a 3D numpy
%array, which was then converted to a .mat file
%import these files
Ann = importdata('Ann.mat');
CT = importdata('CT.mat');
Dose = importdata('Dose.mat');

%plot(Ann(26,:,123)/max(Ann(26,:,123)))
%hold on;
plot(Dose(26,:,123)/max(Dose(26,:,123)))
hold off;
%plot(x_axis_smooth1,y_smoothAnn1/max(y_smoothAnn1)*0.95,Color="red")
%hold on;
%plot(x_axis_smooth2,y_smoothAnn2/max(y_smoothAnn2)*0.95,Color="red")
%hold on;

%get rid of the first and last 98 from x, first 54 and last 49 from y-axis, 
% first 29 and last 25 from dose in z-axis, 
trimDose = Dose(99:151,55:155,30:180);
maxDose = max(trimDose,[],"all")
minDose = min(trimDose,[],"all")

%y = 53
%z=21
%x=120

y = 43
z=26
x=122

%annihilations and dose are flipped along x and y compared to CT
Ann_rot = flip(Ann,3);
Ann_rot = flip(Ann_rot,2);
Dose_rot = flip(Dose,3);
Dose_rot = flip(Dose_rot,2);

save("Annihilations_rot.mat","Ann_rot")
save("Dose_rot","Dose_rot")
save("CT_rot", "CT")

%verify the profiles are consistent based on identical profiles from dicom
%viewer
%plot(squeeze(Ann(z,y,:)))
%hold on;
%plot(squeeze(CT(z,y,:)))
%hold on;
%plot(squeeze(Ann_test(z,y,:)))
%hold off;

%print a 2D histogram of dose vs annihilations
histogram2(Dose_rot,Ann_rot,500);
%histogram2(CT,Ann_rot,250);

%import annihilations to dose data from python and plot profiles
Ann_To_Dose = importdata('Ann_To_Dose.mat');
CT_from_Python = importdata('CT_from_Python.mat');
Dose_from_Python = importdata('Dose_from_Python.mat');

smoothAnn = Ann_To_Dose(z,:,x);
Ann1 = smoothAnn(6:34);
Ann2 = smoothAnn(40:101);
Ann_all = smoothAnn(1:101)
smoothAnn1 = smoothdata(Ann1,"movmedian",5);
smoothAnn2 = smoothdata(Ann2,"movmedian",8);
smoothAnn = smoothdata(Ann_all,"movmedian",10);

N = 2;
N_Ann = 3;

x_axis_smooth = uint32(1+3.25):uint32(length(smoothAnn)+3.25)
x_axis_smooth1 = uint32(6+N):uint32(34+(N));
x_axis_smooth2 = uint32(40+N):uint32(101+(N));
x_axis_ann = uint32(N_Ann):uint32(length(Ann_rot(z,:,x))+(N_Ann-1));
x_axis_dose = uint32(1):uint32(length(Dose_from_Python(z,:,x)));

y_smoothAnn = smoothAnn;
y_smoothAnn1 = smoothAnn1;
y_smoothAnn2 = smoothAnn2;
y_rawAnn = Ann_rot(z,:,x)
y_dose = Dose_from_Python(z,:,x)

CT_profile = cast(squeeze(CT_from_Python(z,:,x)),"single");
CT_profile_pos = CT_profile(1:100)-min(CT_profile(1:100));
CT_profile_norm = CT_profile_pos/max(CT_profile_norm);

%plot(squeeze(Ann_To_Dose(z,y,:)))
%hold on;
%plot(smoothAnn)
%hold on;
%plot(x_axis_smooth1,y_smoothAnn1/max(y_smoothAnn1)*0.95,Color="red")
%hold on;
%plot(x_axis_smooth2,y_smoothAnn2/max(y_smoothAnn2)*0.95,Color="red")
%hold on;

normsmoothAnn = y_smoothAnn/max(y_smoothAnn);
normDose = y_dose/max(y_dose);

error = mean(((normsmoothAnn(9:85)-normDose(9:85))./normDose(9:85))*100)

plot(x_axis_smooth,y_smoothAnn/max(y_smoothAnn),Color="black",LineWidth=1)
hold on;
%plot(x_axis_ann,y_rawAnn/max(y_rawAnn),color="red",LineWidth=1)
%hold on;
%plot(CT_profile_norm,Color="blue");
%hold on;
plot(x_axis_dose,y_dose/max(y_dose),color="magenta",LineWidth=1)
hold off;
ylim([-0.001,1.1])
xlim([0,105])
legend("Converted Dose","MC Dose")
xlabel("Distance (arb)")
ylabel("Normalized Values")
title("Converted Dose vs MC Dose")