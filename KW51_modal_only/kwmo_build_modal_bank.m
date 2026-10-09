function B=kwmo_build_modal_bank(originalRoot,outputFile,nModes)
% Assemble the blueprint FE and audit against the retained static bank.
if nargin<3,nModes=50;end
root=kwmo_setup(originalRoot);old=pwd;clean=onCleanup(@()cd(old)); %#ok<NASGU>
cd(fullfile(root,'dependencies','official_fe'));
S=load('sensor_definition.mat');
B=load(fullfile(root,'fe_models','model_15.mat'));B.xgrid=B.xgrid(:);
[ok,ix]=ismember(B.labels(1:5),S.labels_s);assert(all(ok));
assert(isequal(S.labels_d(:),B.labels(6:7).'));
H=[S.Ss(ix,:)*1e6;S.Sd*1e3];
[K,M,T]=model_kw51(1,1,0,1,1,0,0,1,[]);
K=(K+K')/2;M=(M+M')/2;
[phi,omega,checks]=kwmo_solve_modes(K,M,nModes);
B.frequency_Hz=omega(:)/(2*pi);B.output_modes=full(H*T*phi);
B.sensor_modes=full(S.Sa*T*phi);B.sensor_mode_labels=S.labels_a;
names=fieldnames(checks);for j=1:numel(names),B.(names{j})=checks.(names{j});end
B.modal_force0=zeros(numel(B.xgrid),nModes,2);B.modal_force2=B.modal_force0;
Z=K\(T.'*H.');Q=Z.'*T.';newH0=zeros(size(B.H0));newH2=zeros(size(B.H2));
tracks={'A','B'};
for tr=1:2
 P0=xrail_2_P(tracks{tr},B.xgrid,1,zeros(size(B.xgrid)));
 P30=xrail_2_P(tracks{tr},B.xgrid,1,30*ones(size(B.xgrid)));P2=(P30-P0)/900;
 B.modal_force0(:,:,tr)=full(phi.'*T.'*P0).';
 B.modal_force2(:,:,tr)=full(phi.'*T.'*P2).';
 newH0(:,:,tr)=(Q*P0).';newH2(:,:,tr)=(Q*P2).';
end
% Check each output separately: do not mix strain and displacement units.
% 1e-4 is supported by the supplied-model SciPy/Octave audit (see report).
% It is a numerical replay tolerance, not permission to alter FE parameters.
newG=full(Q*S.Sd.');B.static_channel_relative_errors=zeros(7,5);
for tr=1:2
 B.static_channel_relative_errors(:,tr)=relative_columns(newH0(:,:,tr),B.H0(:,:,tr));
 B.static_channel_relative_errors(:,tr+2)=relative_columns(newH2(:,:,tr),B.H2(:,:,tr));
end
B.static_channel_relative_errors(:,5)=relative_columns(newG.',B.G.');
B.static_bank_relative_error=max(B.static_channel_relative_errors(:));
B.static_bank_tolerance=1e-4;
assert(B.static_bank_relative_error<B.static_bank_tolerance, ...
 'Static/modal FE baseline audit failed: max output-wise relative error %.3g (limit %.3g).', ...
 B.static_bank_relative_error,B.static_bank_tolerance);
B.note='Fixed baseline mode shapes; modal frequency/damping variability. Stored H0/H2/G retained exactly.';
kwpe_save_mat(outputFile,B);
fprintf('Modal bank: mass %.3g; eigen %.3g; static max-channel %.3g\n', ...
 B.mass_orthogonality_error,B.eigen_residual,B.static_bank_relative_error);
end
function r=relative_columns(actual,reference)
den=sqrt(sum(abs(reference).^2,1));num=sqrt(sum(abs(actual-reference).^2,1));
assert(all(den>0),'Unexpected zero reference response in static bank.');r=(num./den).';
end
