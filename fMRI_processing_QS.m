
clear all
spmpath=('/local/spm12-20160115');
addpath('/local/spm12-20160115');
spm_jobman('initcfg');
spm('defaults','FMRI');
matlabbatch = {};




%P{1} = {'ABEV'; 'CLJU'; 'COLE'; 'DAAD'; 'ELAS'; 'EMHO'; 'GAGR'; 'HAAB'; 'JAGO'; 'JEPI'; 'LABI'; 'MAPE'; 'MATO'; 'NAPA'; 'NASA'; 'NIUP'; 'RAWA'; 'STNI'; 'VIPA'; 'WEEV'; 'WIJO'; 'WIPE'};  % controls
%P{2} = {'ATSA'; 'CHEH'; 'CHAK'; 'ILIK'; 'HALO'; 'JOPE'; 'MACO'; 'PAEV'; 'SAEL'; 'SAWH'}; % LTLE
%P{3} = {'DALE'; 'DASC'; 'ELWE'; 'EWDE'; 'FRNI'; 'HAUP'; 'JOGR'; 'RETH'; 'TIPI'}; % RTLE
P{1} = {};  % controls
P{2} = {}; % LTLE
P{3} = {'DASC'}; % RTLE
Time= ['01'];	  %Time (01=preop, 02=4months postop, 03=12months postop)


task = {'Face','WordEnc','WordRetr'};
version = ''; %defined below based on the excel sheet
 

group = [''];	%this variable will be defined automatically below. No need to change here.
subgroup = {'Controls'; 'LTLE'; 'RTLE'};
project_dir = ['/mounts/auto/p1602mem'];
main_dir = [fullfile(project_dir,'Data')];
 
for g=1:numel(subgroup)			%for each subgroup
    if strcmp(subgroup{g},'LTLE')	%define the group
    group = ['Patients'];
    elseif strcmp(subgroup{g},'RTLE')
    group = ['Patients'];
    else
    group = [''];
    disp(group)
    end
    
for s=1:numel(P{g})	%for each subject, within each subgroup
mainpath = [fullfile(main_dir,group,subgroup{g},P{g}{s},Time,'/')];

if exist(mainpath)



T1 = dir(char(fullfile(mainpath, 'T1/', strcat(P{g}{s},'_*.nii')))); %T1
T1b = fullfile(T1.folder,T1.name);
T1path = T1.folder;


Database = ['/mounts/auto/p1602mem/Database/P1602_scaninfo.xlsx']; 
range='A:AZ';

%%%% find the fMRI version for each participant %%%% 
    sheet = 3;
    [num,txt,raw] = xlsread(Database,sheet,range); 
    
    fmri_version = raw(strcmp(raw(:,7),Time) & strcmp(raw(:,3),P{g}{s}),9); 
    fmri_version=string(fmri_version{1});
    fmri_version=char(fmri_version);
    
    if strcmp(fmri_version,'A')	
    version = 'versionA';
        elseif strcmp(fmri_version,'B')	
    version = 'versionB';
        elseif strcmp(fmri_version,'C')	
    version = 'versionC';
    else version = 'versionD';
    end 
    

disp('***starting new participant***')
disp(subgroup{g})
disp(P{g}{s})
disp(version) 	
disp(Time)	
	

%%% add summary of processes in .txt file within each participant's folder %%%
fid=fopen(fullfile(mainpath,'process_summary.txt'), 'w');
fprintf(fid, '%s\n%s\n%s\n%s', P{g}{s}, subgroup{g}, Time, version);
fclose(fid);


%%% generate onsets for event-related analysis %%%
out_folder = [fullfile(project_dir,'Onset')];

if ~isempty(dir(fullfile(out_folder, strcat('Onset_', P{g}{s},Time,'_WordEnc.mat'))));
disp('---Onset file already generated---')
else

if strcmp(version,'versionD')
score_file_dir = [fullfile(project_dir,'Scripts')];
else score_file_dir = [fullfile(project_dir,'Scripts','Individual_steps','2.fMRI_performance')];
end 
score_file = [fullfile(score_file_dir, strcat('2.fMRI_performance_', version, '_', Time, '.xlsx'))];
Subject = strcat(P{g}{s},Time);

%%% Words encoding %%% Sheet 3
inFile = score_file;
range='A:P';
[num,txt,raw] = xlsread(inFile,3,range); 
SuccEnc=raw(strcmpi(raw,Subject)==1,4); %find information for the whole column
SuccEnc=SuccEnc(~cellfun('isempty',SuccEnc)); %remove empty cells
SuccEnc=cell2mat(SuccEnc); %change format
SuccEnc = transpose(SuccEnc); %transpose from vertical to horizontal
UnsEnc=raw(strcmpi(raw,Subject)==1,5);
UnsEnc=UnsEnc(~cellfun('isempty',UnsEnc));
UnsEnc=cell2mat(UnsEnc);
UnsEnc = transpose(UnsEnc);
names{1}=['SuccEnc'];
names{2}=['UnsEnc'];
durations{1}=[0];
durations{2}=[0];
onsets{1}=[SuccEnc];
onsets{2}=[UnsEnc];
outFile = [fullfile(out_folder, strcat('Onset_', Subject, '_WordEnc.mat'))];
save(outFile, 'names', 'durations', 'onsets');

%%% Words retrieval %%% Sheet 3
inFile = score_file;
range='A:P';
[num,txt,raw] = xlsread(inFile,3,range);
Hits=raw(strcmpi(raw,Subject)==1,6);
Hits=Hits(~cellfun('isempty',Hits));
Hits=cell2mat(Hits);
Hits=transpose(Hits);
Misses=raw(strcmpi(raw,Subject)==1,7);
Misses=Misses(~cellfun('isempty',Misses));
Misses=cell2mat(Misses);
Misses=transpose(Misses);
Familiar=raw(strcmpi(raw,Subject)==1,8);
Familiar=Familiar(~cellfun('isempty',Familiar));
Familiar=cell2mat(Familiar);
Familiar=transpose(Familiar);
if isempty(Familiar);
 Familiar = [1];
end 
CorrRej=raw(strcmpi(raw,Subject)==1,9);
CorrRej=CorrRej(~cellfun('isempty',CorrRej));
CorrRej=cell2mat(CorrRej);
CorrRej=transpose(CorrRej); 
if isempty(CorrRej);
 CorrRej = [1];
end 
FalseAl=raw(strcmpi(raw,Subject)==1,10);
FalseAl=FalseAl(~cellfun('isempty',FalseAl));
FalseAl=cell2mat(FalseAl);
FalseAl=transpose(FalseAl); 
if isempty(FalseAl);
 FalseAl = [1];
end 
names{1}=['Hits'];
names{2}=['Misses'];
names{3}=['Familiar'];
names{4}=['CorrRej'];
names{5}=['FalseAl'];
durations{1}=[0];
durations{2}=[0];
durations{3}=[0];
durations{4}=[0];
durations{5}=[0];
onsets{1}=[Hits];
onsets{2}=[Misses];
onsets{3}=[Familiar];
onsets{4}=[CorrRej];
onsets{5}=[FalseAl];
outFile= [fullfile(out_folder,strcat('Onset_', Subject, '_WordRetr.mat'))];
save(outFile, 'names', 'durations', 'onsets');

%%% Face %%% Sheet 4
inFile = score_file;
range='A:P';
[num,txt,raw] = xlsread(inFile,4,range); 
SuccEnc=raw(strcmpi(raw,Subject)==1,4);
SuccEnc=SuccEnc(~cellfun('isempty',SuccEnc));
SuccEnc=cell2mat(SuccEnc);
SuccEnc=transpose(SuccEnc);
UnsEnc=raw(strcmpi(raw,Subject)==1,5);
UnsEnc=UnsEnc(~cellfun('isempty',UnsEnc));
UnsEnc=cell2mat(UnsEnc);
UnsEnc=transpose(UnsEnc);
Hits=raw(strcmpi(raw,Subject)==1,6);
Hits=Hits(~cellfun('isempty',Hits));
Hits=cell2mat(Hits);
Hits=transpose(Hits);
Misses=raw(strcmpi(raw,Subject)==1,7);
Misses=Misses(~cellfun('isempty',Misses));
Misses=cell2mat(Misses);
Misses=transpose(Misses);
if isempty(Misses);
 Misses = [1];
end 
Familiar=raw(strcmpi(raw,Subject)==1,8);
Familiar=Familiar(~cellfun('isempty',Familiar));
Familiar=cell2mat(Familiar);
Familiar=transpose(Familiar);
if isempty(Familiar);
 Familiar = [1];
end 
CorrRej=raw(strcmpi(raw,Subject)==1,9);
CorrRej=CorrRej(~cellfun('isempty',CorrRej));
CorrRej=cell2mat(CorrRej);
CorrRej=transpose(CorrRej);
if isempty(CorrRej);
 CorrRej = [1];
end 
FalseAl=raw(strcmpi(raw,Subject)==1,10);
FalseAl=FalseAl(~cellfun('isempty',FalseAl));
FalseAl=cell2mat(FalseAl);
FalseAl=transpose(FalseAl); 
if isempty(FalseAl);
 FalseAl = [1];
end 
 names{1}=['SuccEnc'];
 names{2}=['UnsEnc'];
 names{3}=['Hits'];
 names{4}=['Misses'];
 names{5}=['Familiar'];
 names{6}=['CorrRej'];
 names{7}=['FalseAl'];
 durations{1}=[0];
 durations{2}=[0];
 durations{3}=[0];
 durations{4}=[0];
 durations{5}=[0];
 durations{6}=[0];
 durations{7}=[0];
 onsets{1}=[SuccEnc];
 onsets{2}=[UnsEnc];
 onsets{3}=[Hits];
 onsets{4}=[Misses];
 onsets{5}=[Familiar];
 onsets{6}=[CorrRej];
 onsets{7}=[FalseAl];
outFile= [fullfile(out_folder,strcat('Onset_', Subject, '_Face.mat'))];
save(outFile, 'names', 'durations', 'onsets');
disp('onset generation done')
end 	
	
	
	
	
%%%%%% PREPROCESS MEMORY %%%%%%
disp('----pre-processing starting----')
find_subject_please = strcat('^' , P{g}{s});

%%% Segment T1 %%% 
iT1seg = length(matlabbatch) + 1;
if ~isempty(dir(fullfile(T1path, 'y*.nii')))
disp('---T1 already segmented---')
else

matlabbatch{iT1seg}.spm.spatial.preproc.channel.vols(1) = {T1b};
matlabbatch{iT1seg}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{iT1seg}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{iT1seg}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spmpath, '/tpm/TPM.nii,1')};
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spmpath, '/tpm/TPM.nii,2')};
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spmpath, '/tpm/TPM.nii,3')};
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spmpath, '/tpm/TPM.nii,4')};
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spmpath, '/tpm/TPM.nii,5')};
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spmpath, '/tpm/TPM.nii,6')};
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{iT1seg}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{iT1seg}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{iT1seg}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{iT1seg}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{iT1seg}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{iT1seg}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{iT1seg}.spm.spatial.preproc.warp.write = [0 1];
spm_jobman('run',matlabbatch);
end

%%% Normalise T1 %%% 
iwT1 = length(matlabbatch) + 1;
if ~isempty(dir(fullfile(T1path, 'w*.nii')))
disp('---T1 already normalised---')
else
def_field = dir(char(fullfile(T1path,'y_*')));
def = fullfile(T1path, def_field.name);
matlabbatch{iwT1}.spm.spatial.normalise.write.subj.def = {def};
matlabbatch{iwT1}.spm.spatial.normalise.write.subj.resample = {T1b};
%matlabbatch{iwT1}.spm.spatial.normalise.write.subj.resample = cfg_dep('File Selector (Batch Mode): Selected Files (find_subject_please)', substruct('.','val', '{}',{iT1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{iwT1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 78 76 85];
matlabbatch{iwT1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{iwT1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{iwT1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch); 
clear matlabbatch; 
matlabbatch = {};
end


%%% functional images %%%

for t=1:numel(task); %for each task
taskpath = [fullfile(mainpath, 'Memory', task{t})];
disp(strcat('---',task{t},'---'))
cd(taskpath) 

if exist(fullfile(taskpath,'nii',strcat(task{t},'.nii')))
disp(strcat('---4D file already created---'))
else 
disp(strcat('---converting 3D files to 4D---'))

%select 3D files
S3D = length(matlabbatch) + 1;
matlabbatch{S3D}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {fullfile(taskpath,'/nii')};
matlabbatch{S3D}.cfg_basicio.file_dir.file_ops.file_fplist.filter = find_subject_please;
matlabbatch{S3D}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%%% convert 3D to 4D %%%
convert = length(matlabbatch) + 1;
matlabbatch{convert}.spm.util.cat.vols = cfg_dep('File Selector (Batch Mode): Selected Files (find_subject_please)', substruct('.','val', '{}',{S3D}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));                                
matlabbatch{convert}.spm.util.cat.name = strcat(task{t},'.nii');
matlabbatch{convert}.spm.util.cat.dtype = 4;
spm_jobman('run',matlabbatch);
clear matlabbatch; 
matlabbatch = {};

%%%move 3D nii files to backup folder %%%
disp('----moving 3D files to backup folder----')
source = fullfile(main_dir, group, subgroup{g}, P{g}{s},Time,'Memory', task{t},'nii');
destination = fullfile(main_dir, group, subgroup{g}, P{g}{s},Time,'Memory', task{t},'nii','nii_backup');
cd(source) 
movefile(strcat('*',P{g}{s},'*'), destination);
end  

%%% select 4D files %%%
S4D = length(matlabbatch) + 1;
matlabbatch{S4D}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {fullfile(taskpath,'/nii')};
matlabbatch{S4D}.cfg_basicio.file_dir.file_ops.file_fplist.filter = strcat(task{t},'.nii');
matlabbatch{S4D}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

spm_jobman('run',matlabbatch);
clear matlabbatch; 
matlabbatch = {};


%%% unwarp (prefix is ud) %%%
if exist(fullfile(taskpath, '/nii', strcat('ud',task{t},'.nii')));
disp('file already unwarped')
else

if exist(fullfile(mainpath, 'DTI'));
disp(strcat('---find DTI scan info for', {' '}, subgroup{g}, {' '}, P{g}{s}, {' '}, 'at time', {' '}, Time, '---'))
U = char(strcat(task{t}, {' '}, 'unwarped with DTI data'));
fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
fprintf(fid, '\n%s', U);
fclose(fid);	

if strcmp(subgroup{g},'Controls')
 if strcmp(Time,'02')
   DTI_01= fullfile(main_dir,group,subgroup{g},P{g}{s},'01', 'DTI');
   if exist(DTI_01)
     DTI_02= fullfile(main_dir,group,subgroup{g},P{g}{s},'02', 'DTI');
     mkdir DTI_02
     copyfile(DTI_01,DTI_02);
     fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
     fprintf(fid, '\n%s', 'DTI folder from time 01 copied to time 02');
     fclose(fid);
   end     	
 end
end 
 
if strcmp(group, 'Patients'); 
    sheet = 1;
    [num,txt,raw] = xlsread(Database,sheet,range); 
    
    if strcmp(Time,'02') | strcmp(Time,'03') | strcmp(Time,'04'); %post-op patients
    exam = raw(strcmp(raw(:,3),P{g}{s}),23);  %DTI session number (column number 31 in Database)
    exam=string(exam{1});
    exam=char(exam);
    NODDI_scan =  raw(strcmp(raw(:,3),P{g}{s}),29); %DTI-rev serie numbet
    NODDI_scan=string(NODDI_scan{1});
    NODDI_scan=char(NODDI_scan);
    revPE_scan = raw(strcmp(raw(:,3),P{g}{s}),30);  %DTI serie number
    revPE_scan=string(revPE_scan{1});
    revPE_scan=char(revPE_scan);
    
    else %pre-op patients     
    sheet = 1;
    [num,txt,raw] = xlsread(Database,sheet,range); 
    exam = raw(strcmp(raw(:,3),P{g}{s}),10);  %DTI session number 
    exam=string(exam{1});
    exam=char(exam);  
    NODDI_scan =  raw(strcmp(raw(:,3),P{g}{s}),17); %DTI-rev serie number
    NODDI_scan=string(NODDI_scan{1});
    NODDI_scan=char(NODDI_scan);
    revPE_scan = raw(strcmp(raw(:,3),P{g}{s}),18);  %DTI serie number
    revPE_scan=string(revPE_scan{1});
    revPE_scan=char(revPE_scan);
    end
    
else %controls
sheet = 2;
[num,txt,raw] = xlsread(Database,sheet,range); 
exam = raw(strcmp(raw(:,3),P{g}{s}),9);  %DTI session number 
exam=string(exam{1});
exam=char(exam);
NODDI_scan =  raw(strcmp(raw(:,3),P{g}{s}),16); %DTI-rev serie number
NODDI_scan=string(NODDI_scan{1});
NODDI_scan=char(NODDI_scan);
revPE_scan = raw(strcmp(raw(:,3),P{g}{s}),17);  %DTI serie number
revPE_scan=string(revPE_scan{1});
revPE_scan=char(revPE_scan);
end 


topup_script= '/home/sbuck/Desktop/p1602mem/Scripts/Other/diffusion_topup_Sjoerd';
[r st] = system([topup_script ' ' exam ' ' NODDI_scan ' ' revPE_scan ' ' taskpath])



elseif exist (fullfile(mainpath, 'Pepolar'))
U = char(strcat(task{t}, {' '}, 'unwarped with Pepolar data'));
fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
fprintf(fid, '\n%s', U);
fclose(fid);	
display('unwarping using pepolar data')
topup_script='/home/sbuck/Desktop/p1602mem/Scripts/Other/distortion_correction_pepolar';
[r st] = system([topup_script ' ' taskpath ' ' P{g}{s}])


else display ('no data found, skipping unwarping')
U = char(strcat(task{t}, {' '}, 'not unwarped'));
fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
fprintf(fid, '\n%s', U);
fclose(fid);	
fullfile(taskpath,'nii',strcat(task{t},'.nii'));
movefile(fullfile(taskpath,'nii',strcat(task{t},'.nii')), fullfile(taskpath,'nii',strcat('ud',task{t},'.nii')));   %add 'ud' prefix to the file anyway, for following processing
end 
end




%%% Realign %%% 
if exist(fullfile(mainpath,'Memory', task{t},'/nii', strcat('meanud',task{t},'.nii')));
disp('file already realigned')
else

ir = length(matlabbatch) + 1;
nii4d = {fullfile(taskpath, '/nii', strcat('ud',task{t},'.nii'))};
matlabbatch{ir}.spm.spatial.realign.estwrite.data = {nii4d};
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{ir}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{ir}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{ir}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{ir}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{ir}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{ir}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
end


%%% Coreg %%% 
if exist(fullfile(mainpath,'Memory', task{t},'/nii', strcat('wrud',task{t},'.nii')));
disp('file already normalised')
else

ic = length(matlabbatch) + 1;
meannii = fullfile(mainpath,'Memory', task{t},'/nii', strcat('meanud',task{t},'.nii')); %meanimage
rnii = fullfile(mainpath,'Memory', task{t},'/nii', strcat('rud',task{t},'.nii')); %realigned images
matlabbatch{ic}.spm.spatial.coreg.estimate.ref = {T1b}; %T1
matlabbatch{ic}.spm.spatial.coreg.estimate.source = {meannii}; %meanimage
matlabbatch{ic}.spm.spatial.coreg.estimate.other = {rnii}; %realigned images
matlabbatch{ic}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{ic}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{ic}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{ic}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];





%%% Normalise %%% 
iw = length(matlabbatch) + 1;
def_field = dir(char(fullfile(T1path,'y_*')));
def = fullfile(T1path, def_field.name);
matlabbatch{iw}.spm.spatial.normalise.write.subj.def = {def};
matlabbatch{iw}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{ic}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{iw}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70 78 76 85];
matlabbatch{iw}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{iw}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{iw}.spm.spatial.normalise.write.woptions.prefix = 'w';
end 




%matlabbatch{iw}.spm.spatial.normalise.estwrite.subj.vol = {T1b}; %T1- this line can be masked out, assuming T1 was already segmented and y_* file already created
%matlabbatch{iw}.spm.spatial.normalise.estwrite.subj.vol = {def};
%matlabbatch{iw}.spm.spatial.normalise.estwrite.subj.resample = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{ic}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'))
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/mounts/auto/linux-local/spm12-20160115/tpm/TPM.nii'};
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
%matlabbatch{iw}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
%matlabbatch{iw}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70 78 76 85];
%matlabbatch{iw}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
%matlabbatch{iw}.spm.spatial.normalise.estwrite.woptions.interp = 4;
%matlabbatch{iw}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';







%%% Smooth %%% 
if exist(fullfile(mainpath,'Memory', task{t},'/nii', strcat('swrud',task{t},'.nii')));
disp('file already smoothed')

else
is = length(matlabbatch) + 1;
wnii = fullfile(mainpath,'Memory', task{t},'/nii', strcat('wrud',task{t},'.nii'));
matlabbatch{is}.spm.spatial.smooth.data = {wnii};
matlabbatch{is}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{is}.spm.spatial.smooth.dtype = 0;
matlabbatch{is}.spm.spatial.smooth.im = 0;
matlabbatch{is}.spm.spatial.smooth.prefix = 's';

spm_jobman('run',matlabbatch);
end

clear matlabbatch; 
matlabbatch = {};














%%%%%% 1ST LEVEL MEMORY %%%%%%
disp('----1st level starting----')

con_table = fullfile('/mounts/auto/p1602mem/Scripts/contrast_files', strcat('contrasts_Block_', task{t}, '.xlsx'));
con_table_event = fullfile('/mounts/auto/p1602mem/Scripts/contrast_files', strcat('contrasts_Event_', task{t}, '.xlsx'));
contrastTable = readtable(con_table);
contrastTable_event = readtable(con_table_event);

if strcmp(task{t},'WordRetr') == 1
onset_block = {fullfile('/home/sbuck/p1602mem/Scripts/Onset_files', strcat('onsets_block_', task{t}, '_', version, '.mat'))};
else
onset_block = {fullfile('/home/sbuck/p1602mem/Scripts/Onset_files', strcat('onsets_block_', task{t}, '.mat'))};
end 

%resection mask for post-op patients and GM mask for the rest%
if strcmp(group, 'Patients');
 if strcmp(Time,'02') | strcmp(Time,'03') | strcmp(Time,'04'); 
  imask = length(matlabbatch) + 1; 
  maskside = regexprep(subgroup{g},'.*/','');
  find_mask = strcat('^' , maskside);
  M = char(strcat(maskside, {' '}, 'resection mask applied for', {' '}, task{t}));
  fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
  fprintf(fid, '\n%s', M);
  fclose(fid);
  disp(M)

  matlabbatch{imask}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/home/sbuck/Desktop/p1602mem/Masks/Resection_masks/'};
  matlabbatch{imask}.cfg_basicio.file_dir.file_ops.file_fplist.filter = find_mask;
  matlabbatch{imask}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
 else 
 iGM = length(matlabbatch) + 1;   
 matlabbatch{iGM}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/home/sbuck/Desktop/p1602mem/Masks/GreyMatter/'};
 matlabbatch{iGM}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^GM';
 matlabbatch{iGM}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
  M = char(strcat('Only grey matter mask applied for', {' '}, task{t}));
  fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
  fprintf(fid, '\n%s', M);
  fclose(fid);
  disp(M)
 end

else
 iGM = length(matlabbatch) + 1;   
 matlabbatch{iGM}.cfg_basicio.file_dir.file_ops.file_fplist.dir = {'/home/sbuck/Desktop/p1602mem/Masks/GreyMatter/'};
 matlabbatch{iGM}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^GM';
 matlabbatch{iGM}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
  M = char(strcat('Only grey matter mask applied for', {' '}, task{t}));
  fid=fopen(fullfile(mainpath,'process_summary.txt'), 'at');
  fprintf(fid, '\n%s', M);
  fclose(fid);
  disp(M)
end


task_folder = fullfile(taskpath, 'nii/');
if exist(task_folder)
task_nii = fullfile(task_folder, strcat(task{t}, '.nii/'));
if exist (task_nii)


%select swru
swru_file = dir(char(fullfile(task_folder, strcat('swrud',task{t},'.nii'))));
swru = fullfile(swru_file.folder,swru_file.name);


%select rp
rp_file = dir(char(fullfile(task_folder, strcat('rp_ud',task{t},'.txt'))));
rp = fullfile(rp_file.folder,rp_file.name);




%%% Block Model Specification %%%
if isdir(fullfile(taskpath,'/Block_unwarp'))
    rmdir(fullfile(taskpath,'/Block_unwarp'), 's'); %if already exists, delete it, to avoid message from SPM "SPM.mat already exists...overwriting etc"
    end
mkdir(fullfile(taskpath,'/Block_unwarp')); %create directory
 

imodel = length(matlabbatch) +1;
matlabbatch{imodel}.spm.stats.fmri_spec.dir = {fullfile(taskpath,'/Block_unwarp')};
matlabbatch{imodel}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{imodel}.spm.stats.fmri_spec.timing.RT = 2.5;
matlabbatch{imodel}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{imodel}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{imodel}.spm.stats.fmri_spec.sess.scans = cellstr(swru);
matlabbatch{imodel}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'pmod', {}, 'orth', {});
matlabbatch{imodel}.spm.stats.fmri_spec.sess.multi = onset_block;
matlabbatch{imodel}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{imodel}.spm.stats.fmri_spec.sess.multi_reg = cellstr(rp);
matlabbatch{imodel}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{imodel}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{imodel}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{imodel}.spm.stats.fmri_spec.volt = 1;
matlabbatch{imodel}.spm.stats.fmri_spec.global = 'None';
matlabbatch{imodel}.spm.stats.fmri_spec.mthresh = 0.8;

if strcmp(group,'Patients'); 
 if strcmp(Time,'02') || strcmp(Time,'03') || strcmp(Time,'04'); 
  matlabbatch{imodel}.spm.stats.fmri_spec.mask = cfg_dep('File Selector (Batch Mode): Selected Files (find_mask)', substruct('.','val', '{}',{imask}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
else matlabbatch{imodel}.spm.stats.fmri_spec.mask = cfg_dep('File Selector (Batch Mode): Selected Files (^GM)', substruct('.','val', '{}',{iGM}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{imodel}.spm.stats.fmri_spec.cvi = 'AR(1)';
 end
else matlabbatch{imodel}.spm.stats.fmri_spec.mask = cfg_dep('File Selector (Batch Mode): Selected Files (^GM)', substruct('.','val', '{}',{iGM}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files')); 
end



%%% Block Model Estimation %%%
iest = length(matlabbatch) + 1;
matlabbatch{iest}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{imodel}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{iest}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{iest}.spm.stats.fmri_est.method.Classical = 1;


%%% Block Contrast Manager %%%
icon = length(matlabbatch) + 1;
matlabbatch{icon}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{iest}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for i = 1:height(contrastTable)
   matlabbatch{icon}.spm.stats.con.consess{i}.tcon.name = contrastTable{i,1}{1};
   matlabbatch{icon}.spm.stats.con.consess{i}.tcon.weights = str2num(contrastTable{i,2}{1});
   %matlabbatch{icon}.spm.stats.con.consess{i}.tcon.convecname = contrastTable{i,2}{1};
   matlabbatch{icon}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
end 

matlabbatch{icon}.spm.stats.con.delete = 1;






%%% Event Model Specification %%%
if isdir(fullfile(taskpath,'/Event_unwarp'))
    rmdir(fullfile(taskpath,'/Event_unwarp'), 's'); %if already exists, delete it, to avoid message from SPM "SPM.mat already exists...overwriting etc"
    end
mkdir(fullfile(taskpath,'/Event_unwarp')); %create directory
 

im_e = length(matlabbatch) + 1;
matlabbatch{im_e}.spm.stats.fmri_spec.dir = {fullfile(taskpath,'/Event_unwarp')};
matlabbatch{im_e}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{im_e}.spm.stats.fmri_spec.timing.RT = 2.5;
matlabbatch{im_e}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{im_e}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
matlabbatch{im_e}.spm.stats.fmri_spec.sess.scans = cellstr(swru);
matlabbatch{im_e}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'pmod', {}, 'orth', {});
matlabbatch{im_e}.spm.stats.fmri_spec.sess.multi = {fullfile('/mounts/auto/p1602mem/Onset/', strcat('Onset_',P{g}{s},Time,'_',task{t},'.mat'))};
matlabbatch{im_e}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{im_e}.spm.stats.fmri_spec.sess.multi_reg = cellstr(rp);
matlabbatch{im_e}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{im_e}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{im_e}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{im_e}.spm.stats.fmri_spec.volt = 1;
matlabbatch{im_e}.spm.stats.fmri_spec.global = 'None';
matlabbatch{im_e}.spm.stats.fmri_spec.mthresh = 0.8;

if strcmp(group,'Patients'); 
 if strcmp(Time,'02') || strcmp(Time,'03') || strcmp(Time,'04'); 
  matlabbatch{im_e}.spm.stats.fmri_spec.mask = cfg_dep('File Selector (Batch Mode): Selected Files (find_mask)', substruct('.','val', '{}',{imask}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
else matlabbatch{im_e}.spm.stats.fmri_spec.mask = cfg_dep('File Selector (Batch Mode): Selected Files (^GM)', substruct('.','val', '{}',{iGM}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{im_e}.spm.stats.fmri_spec.cvi = 'AR(1)';
 end
else matlabbatch{im_e}.spm.stats.fmri_spec.mask = cfg_dep('File Selector (Batch Mode): Selected Files (^GM)', substruct('.','val', '{}',{iGM}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files')); 
end


%%% Event Model Estimation %%%
iest_e = length(matlabbatch) + 1;
matlabbatch{iest_e}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{im_e}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{iest_e}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{iest_e}.spm.stats.fmri_est.method.Classical = 1;


%%% Event Contrast Manager %%%
icon_e = length(matlabbatch) + 1;
matlabbatch{icon_e}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{iest_e}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

for i = 1:height(contrastTable_event)
   matlabbatch{icon_e}.spm.stats.con.consess{i}.tcon.name = contrastTable_event{i,1}{1};
   matlabbatch{icon_e}.spm.stats.con.consess{i}.tcon.weights = str2num(contrastTable_event{i,2}{1});
   matlabbatch{icon_e}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
end 

matlabbatch{icon_e}.spm.stats.con.delete = 1;
end %end "if task folder exists"
end %end "if nii file exists within task folder"

end %end task loop
spm_jobman('run',matlabbatch);
clear matlabbatch; 
matlabbatch = {};
disp(strcat('---1st level memory done for', {' '}, P{g}{s}, '---'))







%%%% move cons %%%% 
disp(strcat('***moving cons for', {' '}, P{g}{s},Time, {' '}, 'to the right folders for 2nd level analyses***'))
analysis = {'Face_Block'; 'Face_Event'; 'WordEnc_Block'; 'WordEnc_Event'; 'WordRetr_Block'; 'WordRetr_Event'};

C{1} = {'Baseline'; 'Face encoding'; 'Face encoding > baseline'; 'Face recognition'; 'Face recognition > baseline'}; % Face block
C{2} = {'Hits'; 'Misses'; 'Correct rejections'; 'Familiar'; 'Hits > Misses'; 'Hits > Familiar'; 'Hits > Misses + Familiar'; 'Hits > Misses + Correct rejections'; 'Hits > Misses + Correct rejections + Familiar'; 'Hits + Correct rejections > Familiar + Misses'; 'Successful > Unsuccessful encoding'; 'Unsuccessful > successful encoding'; 'Successful > unsuccessful'; 'Successful encoding'; 'Unsuccessful encoding'}; % Face event
C{3} = {'Baseline > word encoding'; 'Word encoding'; 'Word encoding > baseline'; 'Baseline'}; % Word Enc block
C{4} = {'Successful encoding'; 'Unsuccessful encoding'; 'Successful > Unsuccessful encoding'; 'Unsuccessful > Successful encoding'}; % Word Enc event
C{5} = {'Baseline'; 'Baseline > Recall'; 'Baseline > Recognition'; 'Recall'; 'Recall > Baseline'; 'Recall > Recognition'; 'Recognition'; 'Recognition > Baseline'; 'Recognition > Recall'}; % Word Retr block
C{6} = {'Hits'; 'Misses'; 'Correct rejections'; 'Familiar'; 'Hits > Misses'; 'Hits > Misses + Correct rejections'; 'Hits > Misses + Correct rejections + Familiar'; 'Hits > Misses + Familiar'; 'Hits > Familiar'; 'Hits > Familiar + Correct rejections'; 'Hits + Correct rejections > Misses + Familiar'; 'Hits > Correct rejections'}; % Word Retr event


for A=1:numel(analysis) %for each analysis
if strcmp(analysis{A},'Face_Block') || strcmp(analysis{A},'WordEnc_Block') || strcmp(analysis{A},'WordRetr_Block');
design = 'Block';
else design = 'Event';
end


for r=1:numel(C{A})
B = extractBefore(analysis{A},"_");
outfolder = fullfile('/home/sbuck/Desktop/p1602mem/Analyses/Cons/',Time, B, design, group, subgroup{g});
if not(exist(fullfile(outfolder,C{A}{r})));
mkdir(fullfile(outfolder,C{A}{r}));
end 

 

contrastTable = (fullfile('/mounts/auto/p1602mem/Scripts/contrast_files', strcat('contrasts_', design, '_', B, '.xlsx')));
sheet = 1;
range='A:C';
[num,txt,raw] = xlsread(contrastTable,sheet,range); 
con_number = raw(strcmp(raw(:,1),C{A}{r}),3); 
con_number=string(con_number{1});
con_number=char(con_number);

if str2num(con_number) > 9
    con = 'con_00';
else con = 'con_000';
end

destination = fullfile('/home/sbuck/Desktop/p1602mem/Analyses/Cons/',Time, B, design, group, subgroup{g}, C{A}{r});
source = fullfile(main_dir, group, subgroup{g}, P{g}{s},Time,'Memory', B, strcat(design,'_unwarp'));
sourcefile = fullfile(source, strcat(con,con_number, '.nii'));
outputfile = fullfile(destination, strcat(P{g}{s}, Time, '_con_', con_number, '.nii'));
if exist (sourcefile)
copyfile(sourcefile, outputfile);
end 

end %for each contrast within each analysis e.g. Face encoding > baseline
end %for each analysis e.g. Face_Block






%%%%%% PREPROCESS LANGUAGE %%%%%%

find_subject_please = strcat('^' , P{g}{s});

%%% AN %%%
AN_folder = fullfile(main_dir, group, subgroup{g}, P{g}{s}, Time,'/Language/AN/nii');
if exist(AN_folder)
AN_nii = dir(char(fullfile(AN_folder, strcat(P{g}{s}, '*0000.nii'))));
if exist(fullfile(AN_nii.folder,AN_nii.name));   
disp('---Auditory Naming---')

AN_swr = dir(char(fullfile(AN_folder, strcat('swr',P{g}{s}, '*0000.nii'))));
if exist(fullfile(AN_swr.folder,AN_swr.name));
disp('AN already done')
else
 
%%% select AN files %%%
iAN = length(matlabbatch) + 1;
matlabbatch{iAN}.cfg_basicio.file_dir.file_ops.file_fplist.dir = cellstr(AN_folder) 
matlabbatch{iAN}.cfg_basicio.file_dir.file_ops.file_fplist.filter = find_subject_please;
matlabbatch{iAN}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';
 
%%% Realign AN %%% 
irAN = length(matlabbatch) + 1;
matlabbatch{irAN}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (find_subject_please)', substruct('.','val', '{}',{iAN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{irAN}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{irAN}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{irAN}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{irAN}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{irAN}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{irAN}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%%% Normalise AN %%% 
iwAN = length(matlabbatch) + 1;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.subj.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{irAN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{irAN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.subj.resample(2) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{irAN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.template = {'/home/sbuck/p1103lang/Andre/Templates/cvNSE60symEPI.img,1'};
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -50 78 76 85];
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{iwAN}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';

%%% Smooth AN %%% 
isAN = length(matlabbatch) + 1;
matlabbatch{isAN}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{iwAN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{isAN}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{isAN}.spm.spatial.smooth.dtype = 0;
matlabbatch{isAN}.spm.spatial.smooth.im = 0;
matlabbatch{isAN}.spm.spatial.smooth.prefix = 's';

%select AN swru
iANswr = length(matlabbatch) + 1;
matlabbatch{iANswr}.cfg_basicio.file_dir.file_ops.file_fplist.dir = cellstr(AN_folder);
matlabbatch{iANswr}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^swr';
matlabbatch{iANswr}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%select AN rp
iANrp = length(matlabbatch) + 1;
matlabbatch{iANrp}.cfg_basicio.file_dir.file_ops.file_fplist.dir = cellstr(AN_folder);
matlabbatch{iANrp}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^rp';
matlabbatch{iANrp}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%%% AN Block Model Specification %%%
iANm = length(matlabbatch) + 1;
AN_block_folder = {fullfile(main_dir, group, subgroup{g}, P{g}{s}, Time,'/Language/AN/Block')};
matlabbatch{iANm}.spm.stats.fmri_spec.dir = AN_block_folder 
matlabbatch{iANm}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{iANm}.spm.stats.fmri_spec.timing.RT = 2.5;
matlabbatch{iANm}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{iANm}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.scans = cfg_dep('File Selector (Batch Mode): Selected Files (^swr)', substruct('.','val', '{}',{iANswr}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(1).name = 'AN';
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(1).onset = [2.25
                                                         62.25
                                                         122.25
                                                         182.25
                                                         242.25];
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(1).duration = 25.5;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(2).name = 'AR';
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(2).onset = [47.25
                                                         90
                                                         167.25
                                                         210
                                                         287.25];
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(2).duration = 12.75;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iANm}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
matlabbatch{iANm}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{iANm}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{iANm}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^rp)', substruct('.','val', '{}',{iANrp}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{iANm}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{iANm}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{iANm}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{iANm}.spm.stats.fmri_spec.volt = 1;
matlabbatch{iANm}.spm.stats.fmri_spec.global = 'None';
matlabbatch{iANm}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{iANm}.spm.stats.fmri_spec.mask = {''};
matlabbatch{iANm}.spm.stats.fmri_spec.cvi = 'AR(1)';


%%% AN Block Model Estimation %%%
iANe = length(matlabbatch) + 1;
matlabbatch{iANe}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{iANm}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{iANe}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{iANe}.spm.stats.fmri_est.method.Classical = 1;


%%% AN Block Contrast Manager %%%
iANcon = length(matlabbatch) + 1;
matlabbatch{iANcon}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{iANe}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{iANcon}.spm.stats.con.consess{1}.tcon.name = 'Auditory naming';
matlabbatch{iANcon}.spm.stats.con.consess{1}.tcon.weights = [1 0];
matlabbatch{iANcon}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{iANcon}.spm.stats.con.consess{2}.tcon.name = 'Auditroy reversed';
matlabbatch{iANcon}.spm.stats.con.consess{2}.tcon.weights = [0 1];
matlabbatch{iANcon}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{iANcon}.spm.stats.con.consess{3}.tcon.name = 'Auditory naming - Auditory reversed';
matlabbatch{iANcon}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
matlabbatch{iANcon}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{iANcon}.spm.stats.con.consess{4}.tcon.name = 'Auditory reversed - Auditory naming';
matlabbatch{iANcon}.spm.stats.con.consess{4}.tcon.weights = [-1 1];
matlabbatch{iANcon}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{iANcon}.spm.stats.con.consess{5}.tcon.name = 'Deac AN';
matlabbatch{iANcon}.spm.stats.con.consess{5}.tcon.weights = -1;
matlabbatch{iANcon}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{iANcon}.spm.stats.con.consess{6}.tcon.name = 'Deact AR';
matlabbatch{iANcon}.spm.stats.con.consess{6}.tcon.weights = [0 -1];
matlabbatch{iANcon}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{iANcon}.spm.stats.con.delete = 0;

spm_jobman('run',matlabbatch);
clear matlabbatch;
matlabbatch = {};
end
end
end


%%% PN %%%

PN_folder = fullfile(main_dir, group, subgroup{g}, P{g}{s}, Time,'/Language/PN/nii');
if exist(PN_folder)
PN_nii = dir(char(fullfile(PN_folder, strcat(P{g}{s}, '*0000.nii'))));
if exist(fullfile(PN_nii.folder,PN_nii.name));   
disp('---Picture Naming---')

PN_swr = dir(char(fullfile(PN_folder, strcat('swr',P{g}{s}, '*0000.nii'))));
if exist(fullfile(PN_swr.folder,PN_swr.name));
disp('PN already done')
else

%%% select PN files %%%
iPN = length(matlabbatch) + 1;
matlabbatch{iPN}.cfg_basicio.file_dir.file_ops.file_fplist.dir = cellstr(PN_folder); 
matlabbatch{iPN}.cfg_basicio.file_dir.file_ops.file_fplist.filter = find_subject_please;
matlabbatch{iPN}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%%% Realign PN %%% 
irPN = length(matlabbatch) + 1;
matlabbatch{irPN}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('File Selector (Batch Mode): Selected Files (find_subject_please)', substruct('.','val', '{}',{iPN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{irPN}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{irPN}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{irPN}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{irPN}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{irPN}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{irPN}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%%% Normalise PN %%% 
iwPN = length(matlabbatch) + 1;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.subj.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{irPN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{irPN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.subj.resample(2) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{irPN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.template = {'/home/sbuck/p1103lang/Andre/Templates/cvNSE60symEPI.img,1'};
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -50 78 76 85];
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{iwPN}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';

%%% Smooth PN %%% 
isPN = length(matlabbatch) + 1;
matlabbatch{isPN}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{iwPN}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{isPN}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{isPN}.spm.spatial.smooth.dtype = 0;
matlabbatch{isPN}.spm.spatial.smooth.im = 0;
matlabbatch{isPN}.spm.spatial.smooth.prefix = 's';


%select PN swr
iPNswr = length(matlabbatch) + 1;
matlabbatch{iPNswr}.cfg_basicio.file_dir.file_ops.file_fplist.dir = cellstr(PN_folder);
matlabbatch{iPNswr}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^swr';
matlabbatch{iPNswr}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%select PN rp
iPNrp = length(matlabbatch) + 1;
matlabbatch{iPNrp}.cfg_basicio.file_dir.file_ops.file_fplist.dir = cellstr(PN_folder);
matlabbatch{iPNrp}.cfg_basicio.file_dir.file_ops.file_fplist.filter = '^rp';
matlabbatch{iPNrp}.cfg_basicio.file_dir.file_ops.file_fplist.rec = 'FPList';

%%% PN Block Model Specification %%%
iPNm = length(matlabbatch) + 1;
PN_block_folder = {fullfile(main_dir, group, subgroup{g}, P{g}{s}, Time,'/Language/PN/Block')};
matlabbatch{iPNm}.spm.stats.fmri_spec.dir = PN_block_folder; 

matlabbatch{iPNm}.spm.stats.fmri_spec.dir = PN_block_folder;
matlabbatch{iPNm}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{iPNm}.spm.stats.fmri_spec.timing.RT = 2.5;
matlabbatch{iPNm}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{iPNm}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.scans = cfg_dep('File Selector (Batch Mode): Selected Files (^swr)', substruct('.','val', '{}',{iPNswr}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(1).name = 'PN';
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(1).onset = [3
                                                         78
                                                         153
                                                         228
                                                         303];
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(1).duration = 27;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(2).name = 'ScPics Count';
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(2).onset = [51
                                                         120
                                                         201
                                                         270
                                                         351];
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(2).duration = 12;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(3).name = 'CFace Count';
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(3).onset = [63
                                                         108
                                                         213
                                                         258
                                                         363];
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(3).duration = 12;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep('File Selector (Batch Mode): Selected Files (^rp)', substruct('.','val', '{}',{iPNrp}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{iPNm}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{iPNm}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{iPNm}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{iPNm}.spm.stats.fmri_spec.volt = 1;
matlabbatch{iPNm}.spm.stats.fmri_spec.global = 'None';
matlabbatch{iPNm}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{iPNm}.spm.stats.fmri_spec.mask = {''};
matlabbatch{iPNm}.spm.stats.fmri_spec.cvi = 'AR(1)';

%%% PN Block Model Estimation %%%
iPNe = length(matlabbatch) + 1;
matlabbatch{iPNe}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{iPNm}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{iPNe}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{iPNe}.spm.stats.fmri_est.method.Classical = 1;

%%% PN Block Contrast Manager %%%
iPNcon = length(matlabbatch) + 1;
matlabbatch{iPNcon}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{iPNe}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{iPNcon}.spm.stats.con.consess{1}.tcon.name = 'Picture Naming';
matlabbatch{iPNcon}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{iPNcon}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{2}.tcon.name = 'Scrambled picture count';
matlabbatch{iPNcon}.spm.stats.con.consess{2}.tcon.weights = [0 1];
matlabbatch{iPNcon}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{3}.tcon.name = 'Cartoon faces count';
matlabbatch{iPNcon}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
matlabbatch{iPNcon}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{4}.tcon.name = 'Picture naming - Scrambled pictures count';
matlabbatch{iPNcon}.spm.stats.con.consess{4}.tcon.weights = [1 -1];
matlabbatch{iPNcon}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{5}.tcon.name = 'Picture naming - Cartoon faces count';
matlabbatch{iPNcon}.spm.stats.con.consess{5}.tcon.weights = [1 0 -1];
matlabbatch{iPNcon}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{6}.tcon.name = 'Picture naming - Scrambled pictures & Cartoon faces count';
matlabbatch{iPNcon}.spm.stats.con.consess{6}.tcon.weights = [2 -1 -1];
matlabbatch{iPNcon}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{7}.tcon.name = 'Picture naming Deact';
matlabbatch{iPNcon}.spm.stats.con.consess{7}.tcon.weights = -1;
matlabbatch{iPNcon}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{8}.tcon.name = 'Scrambled picture count Deact';
matlabbatch{iPNcon}.spm.stats.con.consess{8}.tcon.weights = [0 -1];
matlabbatch{iPNcon}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{9}.tcon.name = 'Cartoon faces count Deact';
matlabbatch{iPNcon}.spm.stats.con.consess{9}.tcon.weights = [0 0 -1];
matlabbatch{iPNcon}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{10}.tcon.name = 'Prog Act to PN';
matlabbatch{iPNcon}.spm.stats.con.consess{10}.tcon.weights = [3 2 1];
matlabbatch{iPNcon}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{11}.tcon.name = 'Prog Act to SPc count';
matlabbatch{iPNcon}.spm.stats.con.consess{11}.tcon.weights = [1 3 2];
matlabbatch{iPNcon}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{12}.tcon.name = 'Prog Deact to PN';
matlabbatch{iPNcon}.spm.stats.con.consess{12}.tcon.weights = [-3 -2 -1];
matlabbatch{iPNcon}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.consess{13}.tcon.name = 'Prog Deact to SPc count';
matlabbatch{iPNcon}.spm.stats.con.consess{13}.tcon.weights = [-1 -3 -2];
matlabbatch{iPNcon}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
matlabbatch{iPNcon}.spm.stats.con.delete = 0;

spm_jobman('run',matlabbatch);
end
end
end

disp(strcat('***done processing', {' '}, P{g}{s},'***'))


end %end if mainpath exist 
end %end subject loop
end %end subgroup loop 






