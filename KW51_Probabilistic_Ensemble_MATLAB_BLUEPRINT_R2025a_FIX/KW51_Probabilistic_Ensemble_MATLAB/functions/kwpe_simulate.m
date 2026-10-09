function [nom,delta,draws]=kwpe_simulate(B,E,F,c,varargin)
%KWPE_SIMULATE Quasi-static ensemble and paired discrepancy under blueprint BC.
%
% B.G must contain exactly two support-force response columns corresponding
% to the longitudinal Leuven A/C bearing coordinates. Unlike the legacy
% four-bearing floating formulation, no zero-net-force projection is used.

t=E.time_s(:);
kin=F.kinematics(:).';
H=kw51_basis(B,E,t,kin);
T=numel(t);O=7;
N=c.candidate_count;M=c.paired_count;nAxles=size(H,3);
G=full(B.G);
nSupport=size(G,2);
assert(nSupport==2, ...
    'Blueprint BC requires exactly two Leuven longitudinal support coordinates.');

states=c.stiffness_states(:);
assert(mod(N,numel(states))==0,'CandidateCount must be divisible by the number of stiffness states.');
kappa=repelem(states,N/numel(states));

if isempty(varargin)
    loads=kwpe_load_draws(N,nAxles,c);
else
    provided=varargin{1};
    loads=provided.loads_100kN;
    kappa=provided.stiffness_states(:);
end
assert(isequal(size(loads),[N nAxles])&&numel(kappa)==N,'Invalid candidate draws.');

nom=zeros(N,T,O);
for o=1:O
    h=reshape(H(:,o,:),T,nAxles);
    nom(:,:,o)=(loads*h.')./kappa;
end

if isempty(varargin)
    pairloads=kwpe_load_draws(M,nAxles,c);
else
    pairloads=provided.paired_loads_100kN;
end
nompair=zeros(M,T,O);
for o=1:O
    h=reshape(H(:,o,:),T,nAxles);
    nompair(:,:,o)=(pairloads*h.')/c.pair_nominal_stiffness;
end

psi=kwpe_smooth_basis(t,c.support_frequency_limit_Hz);
if isempty(varargin)
    factors=exp(c.stiffness_log_sd*randn(M,1)+c.common_compliance_log_sd*randn(M,1));
    recovery=exp(c.strain_recovery_log_sd*randn(M,5));
    kins=repmat(kin,M,1);
    kins(:,1)=kins(:,1).*(1+c.speed_relative_sd*randn(M,1));
    kins(:,2)=kins(:,2)+c.entry_sd_s*randn(M,1);
    kins(:,4)=kins(:,4)+c.acceleration_sd_m_s2*randn(M,1);
    forcecoef=randn(M,size(psi,2),nSupport);
else
    factors=provided.paired_stiffness_factors(:);
    recovery=provided.paired_recovery;
    kins=provided.paired_kinematics;
    forcecoef=provided.support_force_coefficients;
end

assert(isequal(size(pairloads),[M nAxles]),'Invalid paired load draws.');
assert(size(forcecoef,1)==M&&size(forcecoef,2)==size(psi,2)&&size(forcecoef,3)==nSupport, ...
    ['Saved discrepancy draws and FE bank are incompatible. Blueprint runs ' ...
     'must contain exactly two support-force coefficient histories.']);

delta=zeros(M,T,O);
for m=1:M
    forces=psi*reshape(forcecoef(m,:,:),size(psi,2),nSupport)*c.support_force_sd_N;
    % Do NOT subtract the mean across supports. Under blueprint BC the
    % Brussels end anchors longitudinal rigid motion, so common Leuven force
    % content is physically admissible and does not create a free rigid mode.
    perturbed=kw51_basis(B,E,t,kins(m,:));
    yp=reshape(reshape(perturbed,T*O,nAxles)*pairloads(m,:).',T,O);
    scale=c.pair_nominal_stiffness*factors(m);
    yp=yp/scale+(forces*G.')/scale;
    yp(:,1:5)=yp(:,1:5).*recovery(m,:);
    delta(m,:,:)=reshape(yp-reshape(nompair(m,:,:),T,O),1,T,O);
end

draws=struct('loads_100kN',loads,'stiffness_states',kappa, ...
    'paired_loads_100kN',pairloads,'paired_stiffness_factors',factors, ...
    'paired_recovery',recovery,'paired_kinematics',kins, ...
    'support_force_coefficients',forcecoef);
end
