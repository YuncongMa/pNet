function ParcellationInit(LeftCell, RightCell, SubjectsIDs, Repeat_Index, prepDataFile, K, InitializeFolder)

mkdir([InitializeFolder '/Input']);
sbjListFile = [InitializeFolder '/Input/sbjListFile_' num2str(Repeat_Index) '.txt'];
system(['rm ' sbjListFile]);
for i = 1:length(SubjectsIDs)
  cmd = ['echo ' LeftCell{SubjectsIDs(i)} ' >> ' sbjListFile];
  system(cmd);
  cmd = ['echo ' RightCell{SubjectsIDs(i)} ' >> ' sbjListFile];
  system(cmd);
end

SubjectsFolder = '/share/apps/freesurfer/6.0.0/subjects/fsaverage5';
% for surface data
surfL = [SubjectsFolder '/surf/lh.pial'];
surfR = [SubjectsFolder '/surf/rh.pial'];
surfML = [SubjectsFolder '/label/lh.Medial_wall.label'];
surfMR = [SubjectsFolder '/label/rh.Medial_wall.label'];

spaR = 1;
vxI = 1;
ard = 0;
iterNum = 2000;
tNum = 569;
alpha = 1;
beta = 10;
resId = 'Initialize';

outDir = [InitializeFolder '/InitializeRes_' num2str(Repeat_Index)];
deployFuncInit_surf_fs(sbjListFile,surfL,surfR,surfML,surfMR,prepDataFile,outDir,spaR,vxI,ard,iterNum,K,tNum,alpha,beta,resId)

