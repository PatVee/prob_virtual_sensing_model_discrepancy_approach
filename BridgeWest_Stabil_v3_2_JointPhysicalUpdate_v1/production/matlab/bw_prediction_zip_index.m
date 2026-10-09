function A=bw_prediction_zip_index(dataRoot,cacheRoot)
%BW_PREDICTION_ZIP_INDEX Index only the raw LT-W traffic archives + metadata.
% ST-W shaker archives are not required because the 84-state shaker ensemble
% is already frozen in the package.
expected=[{'New folder (2).zip'},arrayfun(@(i)sprintf('ltw24_%d.zip',i),1:4,'UniformOutput',false), ...
    arrayfun(@(i)sprintf('ltw21_%d.zip',i),1:2,'UniformOutput',false)];
assert(exist(dataRoot,'dir')==7,'Data folder does not exist: %s',dataRoot);
if ~exist(cacheRoot,'dir'),mkdir(cacheRoot);end
jr=java.io.File(dataRoot);jc=java.io.File(cacheRoot);A.root=char(jr.getCanonicalPath());A.cache=char(jc.getCanonicalPath());
A.members=struct('archive',{},'member',{},'crc',{},'bytes',{},'name',{});A.archives=struct('name',{},'path',{},'bytes',{});
files=dir(fullfile(A.root,'**','*.zip'));
for i=1:numel(expected)
    hits=find(strcmpi({files.name},expected{i}));
    assert(numel(hits)==1,'Expected exactly one %s below %s; found %d.',expected{i},dataRoot,numel(hits));
    f=files(hits);path=fullfile(f.folder,f.name);A.archives(end+1)=struct('name',f.name,'path',path,'bytes',f.bytes); %#ok<AGROW>
    z=java.util.zip.ZipFile(path);cleanup=onCleanup(@()z.close());en=z.entries();
    while en.hasMoreElements()
        e=en.nextElement();if e.isDirectory(),continue;end
        member=char(e.getName());[~,base,ext]=fileparts(strrep(member,'\','/'));
        if ~any(strcmpi(ext,{'.h5','.csv'})),continue;end
        name=[base ext];row=struct('archive',path,'member',member,'crc',double(e.getCrc()),'bytes',double(e.getSize()),'name',name);
        if isempty(A.members),hit=[];else,hit=find(strcmp({A.members.name},name),1);end
        if isempty(hit),A.members(end+1)=row; %#ok<AGROW>
        else,old=A.members(hit);assert(old.crc==row.crc&&old.bytes==row.bytes,'Conflicting copies of %s.',name);end
    end
    clear cleanup
end
names=string({A.members.name});A.counts=struct('ltw2022',sum(startsWith(names,'LTW_ICETD_')&contains(names,'_2022_')), ...
    'ltw2025',sum(startsWith(names,'LTW_ICETD_')&contains(names,'_2025_')),'pls',sum(startsWith(names,'PLSW_')));
A.FE_policy='Detailed Stabil model + frozen shaker-informed modal discrepancy; no ST-W reread required for passage prediction.';
end
