function Directory_List=GUI_Get_Directory(Base_Directory,Window_Title,Multi_Selection)
% By Yuncong Ma, 10/13/2022
% Select one or multiple folders
% Directory_List=GUI_Get_Directory(Base_Directory,Window_Title,Multi_Selection)
% Directory_List is a cell string
% 
% Repackaged from uigetdir2.m from https://www.mathworks.com/matlabcentral/fileexchange/32555-uigetfile_n_dir-select-multiple-files-and-directories


import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

if nargin == 0 || isempty(Base_Directory) % Allow a null argument.
    Base_Directory = pwd;
end
jchooser = javaObjectEDT('javax.swing.JFileChooser', Base_Directory);
jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(Window_Title);
end
jchooser.setMultiSelectionEnabled(Multi_Selection);

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    if Multi_Selection
        jFile = jchooser.getSelectedFiles();
        Directory_List{size(jFile, 1)}=[];
        for i=1:size(jFile, 1)
    		Directory_List{i} = char(jFile(i).getAbsolutePath);
    	end
    else
        jFile = jchooser.getSelectedFile();
        Directory_List=char(jFile(1).getAbsolutePath);
    end
    
	
elseif status == JFileChooser.CANCEL_OPTION
    Directory_List = [];
else
    error('Error occured while picking file.');
end