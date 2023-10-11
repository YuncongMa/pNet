function out = Load_Mathematica_Mat( filename )
% Mathematica will only save 2D matrices in MAT files
Data=load(filename);

if isstruct(Data)
    Field=fieldnames(Data);
    
    if length(Field)==1
        out=Data.(Field{1});
    else
        out=zeros([size(Data.(Field{1})),length(Field)]);
        for i=1:length(Field)
            out(:,:,i)=Data.(Field{i});
        end
    end
else
    out=Data;
end


end

