function [ YF, FU, FC ] = pawn_cdfs(Y,YY,varargin)
%
% Estimate the unconditional output CDF (i.e. when all inputs vary)
% and the conditional CDFs (when one or more inputs are fixed).
%
% Usage:
% [ YF, FU, FC ] = pawn_cdfs(Y,YY)
% [ YF, FU, FC ] = pawn_cdfs(Y,YY,idx)
%
% Input:
%   Y = output sample for estimating the unconditional CDF   - matrix (N,P)
%  YY = output samples for the conditional CDFs          - cell array (M,n) 
%       The element Y{i,k} is the k-th subsample 
%       of output evaluations obtained by fixing the i-th input 
%       (or group of inputs) to its k-th conditioning value,
%       and it is a matrix of size (NC,P).
% idx = when dealing with multiple output (P>1),                   - scalar
%       this function will actually estimate the conditional and
%       unconditional CDFs of one output only. 'idx' thus is the
%       index of the output to be analyzed (i.e. the column of
%       matrix Y and YY{i,k})                            (default value: 1)
%
% Output:
%  YF = values of y at which the CDFs are given              - vector (P,1)
%  FU = values of the empirical unconditional CDF            - vector (P,1)
%       estimated from Y
%  FC = cell array whose element (i,k) is a              - cell array (M,n)
%       vector (P,1) of values of the empirical 
%       conditional CDF estimated from subsample YY{i,k}
%
% Example:
%
% NU = 150 ;
% NC = 100 ;
% n  = 10  ;
% M  = 3   ;
% XU = AAT_sampling('lhs',M,'unif',[-pi,pi],NU);
% YU = model_evaluation('ishigami_homma_function',XU) ;
% XX = pawn_sampling('lhs',M,'unif',[-pi,pi],n,NC);
% YY = pawn_model_evaluation('ishigami_homma_function',XX) ;
% [ YF, FU, FC  ] = pawn_cdfs(YU,YY) ;

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

[NU,P] = size(Y) ;
[M,n]  = size(YY)  ;
[tmp,P2]=size(YY{1});
if P2~=P; error('''Y'' and ''YY{1}''must have the same number of colums'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set optional arguments to their default values:
idx=1;

% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        idx = varargin{1} ;
        if ~isscalar(idx) ; error('input ''idx'' must be scalar'); end
        if abs(idx-round(idx)); error('''idx'' must be integer'); end
        if idx<1  ; error('input ''idx'' must be a positive integer'); end
        if idx>P2 ; error('''idx'' exceeds the number of columns in Y'); end
    end
end

Y=Y(:,idx);

%%%%%%%%%%%%%%%
% Estimate CDFS
%%%%%%%%%%%%%%%

FC     = cell(M,n) ;
Ymin = min(Y) ;
Ymax = max(Y) ;
for i=1:M; 
    for k=1:n ; 
        Ymin = min( Ymin, min(YY{i,k}) ) ; 
        Ymax = max( Ymax, max(YY{i,k}) ) ; 
    end; 
end
    
deltaY = (Ymax-Ymin)/4000 ;
YF =  [ (Ymin-std(Y,'omitnan')/10) : deltaY : (Ymax+std(Y,'omitnan')/10) ]' ;

% unconditional CDF:
FU = empiricalcdf(Y,YF); 

for i=1:M % for each input factor
    
    for k=1:n
        YYik=YY{i,k}; 
        % Approximate conditional CDF:
        FCik = empiricalcdf(YYik(:,idx),YF) ;  
        % Recover results:
        FC{i,k}=FCik;
        
    end
    
end
