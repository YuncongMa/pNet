function [FN, TC] = GIGICA(Data, gFN, varargin)
% Yuncong Ma, 2/1/2024
% GIGICA method for obtaining pFNs from individual fMRI Data and
% group-level FNs gFN
% FN, TC] = GIGICA(Data, gFN)
% This code is adapted from GIG_GIGICA.m Copyright (c) Yuhui DU
% Data is [Nt Nv]
% gFN is [Nv K]
% Optional inputs include threshold_eign, maxIter,a EGv, ErChuPai, ftol,
% error, Nembda, and Spatial_Correspondence
% Spatial correspondence constrain is integrated to ensure the one-to-one
% match between gFNs and pFNs


Options.threshold_eign=eps(class(Data)); % the original version uses eps
Options.maxIter=100;
Options.a=0.5;
Options.EGv=0.3745672075;
Options.ErChuPai=2/pi;
Options.ftol=0.02;
Options.error=1.e-5;
Options.Nembda=1;
Options.Spatial_Correspondence=1;
Options=fOption('GIGICA',Options,varargin);
if isempty(Options)
    return;
end


% Connect to the orignal code
FmriMatr = Data; % [Nt Nv]
ICRefMax = gFN'; % [k, Nv]

thr = eps(class(FmriMatr)); %you can change the parameter. such as thr=0.02;
[n, m] = size(FmriMatr);
FmriMat = FmriMatr - repmat(mean(FmriMatr,2),[1,m]);
CovFmri = (FmriMat*FmriMat') / m;

% Get eigen vectors and values for the correlation between node-wise signal
[Esort, dsort] = eig(CovFmri);
dsort=diag(dsort);
% correct sign of Esort
Esort = Esort .* heaviside(mean(Esort));

% Remove negative eigenvalues
filter_inds=find(sign(dsort)==-1);
dsort(filter_inds)=[];
Esort(:,filter_inds)=[];
dsort = abs(dsort);
%dsort = diag(dsort);
[dsort, flipped_inds] = sort(dsort, 'descend');
numpc = sum(dsort > thr);
Esort = Esort(:, flipped_inds);

EsICnum = size(ICRefMax, 1);

Epart=Esort(:,1:numpc);
dpart=dsort(1:numpc);
Lambda_part=diag(dpart);
WhitenMatrix = (sqrtm(Lambda_part)) \ Epart';
Y=WhitenMatrix*FmriMat;

if thr<1e-10&&numpc<n
    for i=1:size(Y,1)
        Y(i,:)=Y(i,:)/std(Y(i,:));
    end
end

Yinv=pinv(Y);

% Normalize gFN to zscore
ICRefMaxN=zeros(EsICnum,m);
ICRefMaxC=ICRefMax - repmat(mean(ICRefMax,2),[1,m]);
for i=1:EsICnum
    ICRefMaxN(i,:)=ICRefMaxC(i,:)/std(ICRefMaxC(i,:));
end

% Negative entropy of gFNs
NegeEva=zeros(EsICnum,1);
for i=1:EsICnum
    NegeEva(i)=nege(ICRefMaxN(i,:));
end


YR = (1/m)*Y*ICRefMaxN';

maxIter=Options.maxIter;

a=Options.a;
b=1-a;
EGv=Options.EGv;
ErChuPai=Options.ErChuPai;
ICOutMax=zeros(EsICnum,m);
Nemda=Options.Nembda;
for ICnum=1:EsICnum
    Initialization_Method='Original';
    switch Initialization_Method
        case 'Original'
            reference=ICRefMaxN(ICnum,:);
            wc=(reference*Yinv)';
            wc=wc/norm(wc);
            
        case 'Double_Inverse'
            reference=ICRefMaxN(ICnum,:);
            wc=pinv(Y*pinv(ICRefMaxN))';
            wc=wc(:,ICnum);
            wc=wc/norm(wc);
    end
    if Options.Spatial_Correspondence==1
        Source=wc'*Y;
        [~, ps]= max(corr(Source', gFN));
        if sum(ps~=ICnum)>0
            fprintf('Warning: the initial pFN %d does not meet spatial correspondence constrain\n', ICnum);
        end
    end
    
    y1=wc'*Y;
    EyrInitial=(1/m)*(y1)*reference';
    NegeInitial=nege(y1);
    c=(tan((EyrInitial*pi)/2))/NegeInitial;
    IniObjValue=a*ErChuPai*atan(c*NegeInitial)+b*EyrInitial;
   
    % Store history of sources and whether it violates spatial
    % correspondence
%     if Options.Spatial_Correspondence==1
%         History_Source=zeros(maxIter, m);
%         History_SC=zeros(maxIter,1);
%     end

    for i=1:maxIter
        Cosy1=cosh(y1);
        logCosy1=log(Cosy1);
        EGy1=mean(logCosy1);
        Negama=EGy1-EGv;
        tanhy1 = tanh(y1);
        EYgy=(1/m)*Y*(tanhy1)';       
        %EYgy= sum(bsxfun(@times, Y,tanh(y1)/m),2);
        Jy1=(EGy1-EGv)^2;
        KwDaoshu=ErChuPai*c*(1/(1+(c*Jy1)^2));
        %Simgrad=(1/m)*Y*reference';
        Simgrad = YR(:,ICnum);
        g=a*KwDaoshu*2*Negama*EYgy+b*Simgrad;
        d=g/norm(g);
        wx=wc+Nemda*d;
        wx=wx/norm(wx);
        y3=wx'*Y;
        PreObjValue=a*ErChuPai*atan(c*nege(y3))+b*(1/m)*y3*reference';
        ObjValueChange=PreObjValue-IniObjValue;
        ftol=Options.ftol;
        dg=g'*d;
        ArmiCondiThr=Nemda*ftol*dg;
        % Spatial correspondence constraint
%         if Options.Spatial_Correspondence==1
%             Source=wx'*Y;
%             [~, ps]= max(corr(Source', gFN));
%             if sum(ps~=ICnum)>0
%                 %display([ICnum, i, ps])
%                 %wx=wc;
%                 %fprintf('Meet QC constraint\n');
%                 %break
%             else
%                 History_SC(i)=i;
%                 History_Source(i,:)=Source;
%             end
%         end
        % Stop when the change of wx is small enough
        if (wc-wx)'*(wc-wx) <Options.error
            break;
        end

        if ObjValueChange<0
            Nemda=Nemda/2;
            continue;
        end
        IniObjValue=PreObjValue;
        y1=y3;
        wc=wx;
        if ObjValueChange<ArmiCondiThr
            Nemda=Nemda/2;
            continue;
        end
    end
    % find the last source in the history that meets the spatial
    % correspondence constrain
%     if Options.Spatial_Correspondence==1
%         if History_SC(i)>0
%             %display([1,i]);
%             Source=History_Source(i,:);
%         else
%             temp=History_SC(History_SC>0);
%             if isempty(temp)
%                 fprintf('Warning: pFN %d is as same as its corresponding gFN\n', ICnum);
%                 Source=ICRefMaxN(ICnum,:);
%             else
%                 ps=max(temp);
%                 fprintf('Warning: pFN %d is chosen from a previous iteration with spatial correspondence constrain\n', ICnum);
%                 Source=History_Source(ps(end),:);
%             end
%         end
%     else
%         Source=wx'*Y;
%     end
    ICOutMax(ICnum,:)=Source;
end
TCMax=(1/m)*FmriMatr*ICOutMax';

% Check final spatial correspondence once more
SC = corr(ICOutMax', gFN);
Delta_SC = diag(SC) - max(SC-diag(diag(SC)),[],2);
if min(Delta_SC<0)
    fprint('Found %d personalized FNs are not matched to their group-level counterparts\n', sum(Delta_SC<0));
end

FN = ICOutMax';
TC = TCMax;
end


function negentropy=nege(x)

y=log(cosh(x));
E1=mean(y);
E2=0.3745672075;
negentropy=(E1- E2)^2;

end