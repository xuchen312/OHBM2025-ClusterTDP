%--------- Compute approximate TDP bounds using TDPClusters

%--- a) Specify CFT
cft = 0.001;
if ~exist('SPM','var') | isempty(SPM)
    load(spm_select(1,'SPM.mat','Select SPM.mat'));
end
if length(SPM.xCon) > 1
    c = spm_input('Select contrast (index):','+0','i',1,1);
else
    c = 1;
end
switch SPM.xCon(c).STAT
    case 'T'
        df = [1 SPM.xX.erdf];
    case 'F'
        df = [rank(SPM.xCon(c).c) SPM.xX.erdf];
    case 'X'
        df = SPM.xX.erdf;
    otherwise
        error('Unsupported statistic type.');
end
cft = spm_u(cft, df, SPM.xCon(c).STAT);  % convert p-value to statistic

%--- b) Extract k-value from ClusterTDP-SPM output (minClusSz - 1)
kval = 129-1;

%--- c) Generate cluster info files
V = spm_vol(fullfile(SPM.swd, SPM.xCon(1).Vspm.fname));
Tmap = spm_read_vols(V);

% apply CFT & generate clusters
Tmap(Tmap<cft) = 0;
[L,NUM] = spm_bwlabel(Tmap);
% iterate through clusters & export cluster info files
for i = 1:NUM
    if sum(L(:)==i) > kval  % extract clusters with size > kval
        % get [X,Y,Z] voxel coordinates for each cluster
        [x,y,z] = ind2sub(size(L),find(L==i));
        % export to CSV
        fileID = fopen(sprintf('clus%d_size%d_kval%d.csv',i,sum(L(:)==i),kval),'w');
        fprintf(fileID, 'x,y,z\n');
        for j = 1:length(x)
            fprintf(fileID,'%d,%d,%d\n',x(j),y(j),z(j));
        end
        fclose(fileID);
    end
end
