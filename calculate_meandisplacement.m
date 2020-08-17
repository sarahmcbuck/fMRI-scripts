clear all 
addpath(genpath('/local/spm12/toolbox/mw_mfp'));
addpath(genpath('/local/spm12-20160115'));


P{1} = {'ABEV'; 'CLJU'; 'COLE'; 'DAAD'; 'ELAS'; 'EMHO'; 'GAGR'; 'HAAB'; 'JAGO'; 'JEPI'; 'LABI'; 'MAPE'; 'MATO'; 'NAPA'; 'NASA'; 'NIUP'; 'RAWA'; 'STNI'; 'VIPA'; 'WEEV'; 'WIJO'; 'WIPE'};  % controls
P{2} = {'ATSA'; 'CHEH'; 'CHAK'; 'ILIK'; 'HALO'; 'JOPE'; 'MACO'; 'PAEV'; 'SAEL'; 'SAWH'}; % LTLE
P{3} = {'DALE'; 'ELWE'; 'EWDE'; 'FRNI'; 'HAUP'; 'JOGR'; 'RETH'; 'TIPI'}; % RTLE
Time= ('01');	  %Time (01=preop, 02=4months postop, 03=12months postop)

group = ('');	%this variable will be defined automatically below.
subgroup = {'Controls'; 'LTLE'; 'RTLE'};
main_folder = ('/home/sbuck/Desktop/p1602mem/Data/');
out_folder = (fullfile('/home/sbuck/Desktop/p1602mem/Analyses/mean_displacement',Time)); %where the output files will be saved
task = {'Face','WordEnc','WordRetr'};

for t=1:numel(task) %for each task
all_mean_disp = [];

for g=1:numel(subgroup)			%for each subgroup
    if strcmp(subgroup{g},'LTLE')	%define the group
    group = ('Patients');
    elseif strcmp(subgroup{g},'RTLE')
    group = ('Patients');
    else
    group = ('');
    disp(group)
    end
    
for s=1:numel(P{g})	%for each subject, within each subgroup
mainpath = (fullfile(main_folder,group,subgroup{g},P{g}{s},Time));


taskpath = (fullfile(mainpath, 'Memory', task{t}));    
p_data = char(fullfile(taskpath,'nii'));
cd(p_data);

%%% calculate motion fingerprint for each subject 
if ~isempty(dir(fullfile(p_data, 'mw_anamot_results.txt'))) %skip if the motion fingerprint was already calculated
    disp('---motion already estimated---')
else
mw_mfp(pwd,1,0,0) %create motion fingerprint graphs for each subj
mw_anamot(p_data) %extrac the motion fingerprint value for each sub
end 



%%% create a vector with every subject's motion fingerprint (one vector per group)
data = fullfile(p_data,'mw_anamot_results.txt');
T = readtable(data);
mean_disp = T.TD_Mean(1);
all_mean_disp = [all_mean_disp,mean_disp];
outfile_disp = fullfile(out_folder, strcat('mean_displacement_', task{t},'.mat'));
save(outfile_disp)
end
end
end




%%% create plots (one per task
for t=1:numel(task) %for each task
cd(out_folder);
all_disp = strcat('mean_displacement_',task{t},'.mat');
load(all_disp);
plot(all_mean_disp,'color','blue') %plot the data, with blue line
xlabel('\bf Subjects'); 
plot_title = strcat('Total mean displacement for',{' '}, task{t});
title(plot_title)
ylabel('\bf Total mean displacement (mm)'); %in bold

P_1=P{1}';
P_2=P{2}';
P_3=P{3}';
subjects = [P_1,P_2,P_3];
num_subj = numel(subjects); %how many subjects in total
range = 1:num_subj;
set(gca,'XTick',range,'XTickLabel',subjects); %add subject's ID as X labels
xtickangle(45); %tilt the subject's ID
hline(3, 'r'); %add a horizontal line at y=3
ylim([0 4])

cd(out_folder);
fig_title = (strcat('mean_displacement_',task{t},'.fig'));
saveas(gcf,fig_title)
disp(strcat('***plot', {' '},'"',plot_title,'"', {' '}, 'is saved in "analyses" folder***'))
end


%%% combine vectors to create table with all displacements, across tasks 
P = subjects';
T = cell2mat(P);

Table = table(T); %conver matrix to table
cd(out_folder);

load('mean_displacement_Face.mat'); %load each .mat file
Face = (all_mean_disp');
value_F = num2cell(Face);
Table(:,2) = value_F;

load('mean_displacement_WordEnc.mat');
WordEnc = (all_mean_disp');
value_WE = num2cell(WordEnc);
Table(:,3) = value_WE;
load('mean_displacement_WordRetr.mat');
WordRetr = (all_mean_disp');
value_WR = num2cell(WordRetr);
Table(:,4) = value_WR;

Table.Properties.VariableNames(1:4) = {'Participants', 'Face','WordEnc','WordRetr'};
num_subj = numel(subjects); %how many subjects in total
range = 1:num_subj;

y1 = Table.Face;
y2 = Table.WordEnc;
y3 = Table.WordRetr;

%%% create subplot (one plot containing the 3 subplots)
for i=1:3
    subplot(3,1,i)
    if i == 1
    y = y1;
    elseif i == 2
        y = y2;
    else y = y3;
    end
            
    plot(range,y)
    set(gca,'XTick',range,'XTickLabel',subjects); %add subject's ID as X labels
    ylabel(task{i});
    set(get(gca,'ylabel'),'FontWeight','bold','rotation',0, 'Position', [-2.5 1.8 0]); %define properties of y label
    xtickangle(45); %tilt the subject's ID
    hline(3, 'r'); %add a horizontal line at y=3
    ylim([0 4]);
    set(gca,'YTick',[0:1:4]); %change interval y-axis

    %add vertical lines
    text(1,3.4,'Controls'); %add name 
    g_2 = numel(P{1}); %where LTLE category starts
    vline(g_2, 'LineStyle', '--', 'Color', [0.25 0.25 0.25]); %add a v line at the end of the controls category
    text(g_2,3.4,'LTLE'); %add name   
    g_3 = g_2+numel(P{2});
    vline(g_3, 'LineStyle', '--', 'Color', [0.25 0.25 0.25]); %add a v line at the end of the LTLE category
    text(g_3,3.4,'RTLE'); %add name   

end 

xlabel('\bf Participants','FontSize', 12); % x-axis label
suptitle('\bf Total mean displacement (mm)'); % Plot title


%%% save graph
cd(out_folder);
fig_title = ('mean_displacement_all.fig');
saveas(gcf,fig_title)
disp('***subplot is saved in "analyses" folder***')

