function [X,Xf,Xnames,U] = generate_spm(names, onsets, durations, pmod, nscans, TR, nbf, volterra)
% Creates a simulated SPM structure without using the hard disk
%
% names         - {1 x n} cell array of condition names for n conditions
% onsets        - {1 x n} cell array of onsets
% durations     - {1 x n} cell array of durations 
% pmod(i).name  - {1 x p} cell array of p parametric mods for condition i
% pmod(i).param - {1 x p} cell array of p parametric regressors
% nscans        - number of scans
% TR            - volume repetition time
% nbf           - number of haemodynamic basis functions (1,2 or 3)
% volterra      - 1 for 1st order model (default), 2 for 2nd order
%
% Peter Zeidman
% _________________________________________________________________________

if nargin < 8
    volterra = 1;
end

% Build dummy input structure (see spm_get_ons.m)
U = struct();
for i = 1:length(names)
    % Names, onsets, durations, orthogonalization 
    U(i).name = names(i);
    U(i).ons  = onsets{i};
    U(i).dur  = durations{i};
    U(i).orth = false;
    
    % Parametric regressor(s)
    if isempty(pmod) || i > length(pmod) || isempty(pmod(i))
        U(i).P    = struct('name','none','h',[],'i',[],'p',[]);
    else
        % Structure in the right format for for spm_get_ons.m
        P = struct();
        for j = 1:length(pmod(i).name)
            P(j).name = pmod(i).name{j};
            P(j).h = pmod(i).poly{j};
            P(j).P = pmod(i).param{j}(:);
        end
        U(i).P = P;        
    end
    
end

% Basis function spec (see spm_get_bf.m)
xBF = struct();
xBF.UNITS    = 'secs';
xBF.T        = 16;
xBF.T0       = 8;
xBF.Volterra = volterra;
xBF.dt       = TR/xBF.T;
switch nbf
    case 1
        xBF.name = 'hrf';     
    case 2
        xBF.name = 'hrf (with time derivative)';
    case 3
        xBF.name = 'hrf (with time and dispersion derivatives)';
    otherwise
        error('Unknown number of basis functions');
end
xBF = spm_get_bf(xBF);

% Build dummy SPM structure (see spm_get_ons.m)
SPM = struct();
SPM.Sess.U = U;
SPM.xBF    = xBF;
SPM.nscan  = nscans;

% Generate stimulus timeseries (inputs)
U = spm_get_ons(SPM,1);

% For parametric designs with 2nd order Volterra kernels, just retain the 
% parameteric regressor and discard the main regressor (as spm_Volterra.m
% only generates interaction terms for one regressor per condition)
for i = 1:length(U)
    if volterra == 2 && ~isempty(U(i).P.P)
        U(i).name = U(i).name(2);
        U(i).u    = U(i).u(:,2);
    end
end

% Remove spikes from U.u (if there are no parametric regressors) and store
for i = 1:length(U)
    if volterra == 2 && ~isempty(U(i).P.P)       
        maxlevel = max(U(i).P.P);
        U(i).u(U(i).u > maxlevel) = maxlevel;
    else
        U(i).u(U(i).u > 1) = 1;
    end
end

SPM.Sess(1).U = U;

% New design matrix (from spm_fMRI_design)
[X,Xnames] = spm_Volterra(U, SPM.xBF.bf, volterra);

% Resample to scan times (from spm_fMRI_design)
fMRI_T  = SPM.xBF.T;  % microtime bins per scan
fMRI_T0 = SPM.xBF.T0; % microtime offset
X = X((0:(nscans - 1))*fMRI_T + fMRI_T0 + 32,:);

% Add intercept
X = [X ones(nscans,1)];
Xnames{end+1} = 'Intercept';

% High pass filter (128s cut-off)
% K = struct('HParam', 128,...
%              'row', 1:nscans,...
%              'RT',  TR);
% Xf = spm_filter(K,X);
Xf = X;