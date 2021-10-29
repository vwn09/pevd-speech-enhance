function y = pevd_enhance(x,fs,pp)
    %% PEVD_ENHANCE performs PEVD-based speech enhancement [Y]=(X,FS,PP)

    % Usage: y = pevd_enhance(x,fs); % enhance speech using default parameters [1]

    % Inputs (dimensions)
    %   x: time-domain speech (N samples, M channels)
    %   fs: sampling frequency
    %   pp: optional parameters, see default values below

    % Output:
    %   y: decorrelated outputs (N' samples, M channels) -> Enhanced signal is in the first channel

    % The algorithm is controlled by the following parameters [default values]
    %   pp.L: number of iterations [500]
    %   pp.T: frame size [0.1*fs]
    %   pp.W: window length [0.1*fs]
    %   pp.mu: trim factor [10^-3]
    %   pp.delta: max. off-diagonal element [sqrt(N1/M)*10^-2]

    % Author:   Vincent W. Neo
    % Date:     26-10-2021
    
    % Copyright (C) 2017-2021 Vincent W. Neo
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <https://www.gnu.org/licenses/>.

    % References
    % [1]   V. W. Neo, C. Evers, and P. A. Naylor
    %       Enhancement of noisy reverberant speech using polynomial matrix eigenvalue decomposition
    %       IEEE/ACM Trans. Audio, Speech and Lang. Process., vol. 28, 2021. doi: 10.1109/TASLP.2021.3120630
    % [2]   S. Redif, S. Weiss, and J. G. McWhirter
    %       Sequential matrix diagonalisation for polynomial EVD of parahermitian matrices.
    %       IEEE Trans. Signal Process., vol. 63, no. 1, pp. 81-89, Jan. 2015. doi: 10.1109/TSP.2014.2367460
    % [3]   J. G. McWhirter, P. D. Baxter, T. Cooper, S. Redif, and J. Foster
    %       An EVD algorithm for para-hermitian polynomial matrices
    %       IEEE Trans. Signal Process., vol. 55, no. 5, pp. 2158-2169, May 2007. doi: 10.1109/TSP.2007.893222

    %% default parameters
    qq.L = 500;
    qq.mu = 10^-3;
    qq.T = 0.1 * fs;
    qq.W = 0.1 * fs;

    if(nargin>2 && ~isempty(pp))
        pname = fieldnames(pp);
        for ii = 1:length(pname)
            if isfield(pp,pname{ii})
                qq.(pname{ii}) = pp.(pname{ii});
            end
        end
    end
    % extract parameters for algorithm
    maxIt = qq.L;
    mu = qq.mu;
    T = qq.T;
    W = qq.W;

    %% initalisation
    if(size(x,2)>size(x,1)) % to check if input signal has been transposed
        warning('Check audio file dimensions. Processing its transpose...');
        x = x';
    end
    [N, M] = size(x);
    Xfilt = zeros(M,1,N);
    for mm = 1:M
        Xfilt(mm,1,:) = x(:,mm);
    end

    R = frame_stcov(x,T,W); % construct space-time covariance
    [N1,N2,N3,N4] = pmeasure(R);

    if(isfield(qq,'delta'))
        delta = qq.delta;
    else
        delta = sqrt(N1/M)/100; % max. off-diagonal norm
    end

    disp('Computing polynomial eigenvalues and eigenvectors using SMD...');
    [U,D] = smd(R,maxIt,delta,mu);

    disp('Filtering through eigenvector filterbank...');
    y = pmdpdt(U,Xfilt);

    disp('Generating strongly decorrelated output...');
    y = transpose(squeeze(y));
    disp('Decorrelated outputs generated. Enhanced signal in first channel.');
    disp('Done.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sequential Matrix Diagonalization Algorithm Reference Implementation of [2]

function [U,D] = smd(R, maxIt, delta, mu)
    %% SMD performs sequential matrix diagonalisation algorithm for PEVD [2]
    % [U,D] = SMD(R,MAXIT,DELTA,MU)

    % Usage: [U,D] = smd(R,maxit,delta,mu)

    % Inputs
    %   R: parahermitian polynomial matrix (M channels, M channels, W lags)
    %   maxIt: maximum number of iterations
    %   delta: maximum off-diagonal euclidean norm
    %   mu: trim factor

    % Outputs
    % U: paraunitary polynomial matrix with eigenvectors on the rows (not columns) of U
    % D: diagonal polynomial matrix with eigenvalues on the diagonals of D

    %% Initialisation
    it = 1;
    [n1,n2,n3,n4] = pmeasure(R);
    Lc = (size(R,3)+1)/2;

    %% Instantaneous decorrelation at z^0 plane
    [Vs,~] = ordered_evd(R(:,:,Lc));
    U = Vs'; % In SBR2, eigenvectors lie on the rows but standard EVD has them in the columns
    for pp = 1:size(R,3)
        R(:,:,pp) = Vs' * R(:,:,pp) * Vs;
    end

    %% start iteration of SMD algorithm
    [g,cc,ll] = search_max_col(R);

    while(g >= delta && it <= maxIt)
        % Step 1: Search
        R = 0.5*(R + parahermitian(R)); % for numerical robustness

        [g,cc,ll] = search_max_col(R); % ll >= Lc
        Lc = (size(R,3)+1)/2;

        % Step 2: Delay
        R = delay_polymat(R,cc,ll);
        U = delay_ospolymat(U,cc,ll-Lc);

        % Step 3: Zero
        Lc = (size(R,3)+1)/2;
        [Vs,Ds] = ordered_evd(R(:,:,Lc));

        for pp = 1:size(R,3)
            R(:,:,pp) = Vs' * R(:,:,pp) * Vs;
        end

        for pp = 1:size(U,3)
            U(:,:,pp) = Vs' * U(:,:,pp);
        end

        % Step 4: Trim
        if(mu>0)
            thres = mu * n4 / 2;
            R = trim_phmat(R, thres);   % one-sided threshold provided
            U = trim_nphmat(U,mu*size(R,1)/2);
        else
            R = trim_zeros_phmat(R);
        end
        it = it + 1; % increment iteration count
    end
    disp(['PEVD converged after ', num2str(it-1), ' iterations.']);
    D = R;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utility functions for SMD

function [g,j,l] = search_max_col(R)
    % search for the column (less diagonal elements) with the highest norm-2
    [I,J,L] = size(R);
    Lc = 0.5 * (L+1);
    R_hat = off_diag(R);
    Rl = squeeze(sum(abs(R_hat).^2,1));
    [g,I] = max(abs(Rl(:))); % (14) in SMD
    g = sqrt(g);
    [j,l] = ind2sub(size(Rl),I);
end

function A = off_diag(R)
    % returns a polynomial matrix with its diagonal elements = 0 for all lags
    D = zeros(size(R));
    for ii = 1:size(R,3)
        D(:,:,ii) = diag(diag(R(:,:,ii)));
    end
    A = R - D;
end

function A = delay_polymat(R,cc,ll)
    % returns a parahermitian polynomial matrix with the cc-th row and column similarly delayed, (11) in [2]
    [I,J,L] = size(R);
    Lc = (L+1)/2; % z^0-plane
    Le = ll - Lc;

   if(Le>0)
        % shift cc-th column towards z^0 => reduction in polynomial order
        A = cat(3,zeros(I,J,Le),R);
        A(:,cc,1:L) = A(:,cc,Le+1:end);
        A(:,cc,L+1:end) = 0;

        % shift cc-th row towards z^0 => increase in polynomial order
        A = cat(3,A,zeros(I,J,Le));
        A(cc,:,Le+1:end) = A(cc,:,1:end-Le);
        A(cc,:,1:Le) = 0;
    else
        % shift cc-th column backwards => increase in polynomial order
        A = cat(3,R,zeros(I,J,-Le));
        A(:,cc,1-Le:end) = A(:,cc,1:end+Le);
        A(:,cc,1:-Le) = 0;

        % shift cc-th row forward => decrease polynomial order
        A = cat(3,zeros(I,J,-Le),A);
        A(cc,:,1:end+Le) = A(cc,:,1-Le:end);
        A(cc,:,end+Le+1:end) = 0;
    end
end

function A = delay_ospolymat(R,cc,ll)
    % returns a polynomial matrix U which has its cc-th row time advanced/delayed by ll
    [I,J,L] = size(R);
    if(ll>0)
        % want to increase the polynomial order by Le and shift cc-th row
        A = cat(3,R,zeros(I,J,ll));
        A(cc,:,ll+1:end) = A(cc,:,1:L);
        A(cc,:,1:ll) = 0;
    else
        % want to decrease the polynomial order by Le and shift cc-th row
        A = cat(3,zeros(I,J,-ll),R);
        A(cc,:,1:end+ll) = A(cc,:,1-ll:end);
        A(cc,:,end+ll+1:end) = 0;
    end
end

function [Vs,Ds] = ordered_evd(R)
    % performs an ordered EVD notation with eigenvectors lying on the columns of Vs
    [V,D] = eig(R);
    [d,ind] = sort(diag(D),'descend');
    Ds = D(ind,ind);
    Vs = V(:,ind);
end

function R_trim = trim_phmat(R, thres)
    % symmetrically truncate terms in a parahermitian polynomial matrix R
    % if the accumulated F-norm.^2 is less than thres
    ii = 0;
    engy = 0;
    while(engy<thres)
        engy = engy + norm(R(:,:,ii+1),'fro').^2;
        ii = ii + 1;
    end
    R_trim = R(:,:,ii:end-ii+1);
end

function A_trim = trim_nphmat(A, thres)
    % truncate terms in a polynomial matrix U if the accumulated F-norm.^2
    % is less than thres
    [p,q,l] = size(A);
    Lc = (l+1)/2;
    ii = 0;
    jj = 0;
    engy = 0;
    % truncate front
    while(engy<thres && ii<Lc)
        engy = engy + norm(A(:,:,ii+1),'fro').^2;
        ii = ii + 1;
    end
    % truncate back
    engy = 0;
    while(engy<thres && jj<Lc)
        engy = engy + norm(A(:,:,end-jj),'fro').^2;
        jj = jj + 1;
    end
    A_trim = A(:,:,ii:end-jj+1);
end

function R_trim = trim_zeros_phmat(R)
    % if trim=0, we remove planes of zeros on both ends
    ii = 0;
    engy = 0;
    while(engy==0)
        engy = engy + norm(R(:,:,ii+1),'fro');
        ii = ii + 1;
    end
    R_trim = R(:,:,ii:end-ii+1);
end



%% pmatrix utility functions

function [n1,n2,n3,n4] = pmeasure(R)
    % polynomial matrix measures N1-N4 defined in (2.18) in [3]
    n4 = 0;
    for ii = 1:size(R,3)
        n4 = n4 + norm(R(:,:,ii),'fro').^2;
    end
    Lc = (size(R,3)+1)/2; % z^0-plane
    n1 = norm(diag(R(:,:,Lc)),2).^2;
    n2 = norm(R(:,:,Lc),'fro').^2;
    n3 = n2-n1;
end

function C = pmdpdt(A,B)
    % product of 2 polynomial matrices
    [p1,q1,l1] = size(A);
    [q2,r2,l2] = size(B);
    C = zeros(p1,r2,l1+l2-1);
    if(q1==q2)
        for ii = 1:p1
            for jj = 1:r2
                for kk=1:q1
                    C(ii,jj,:) = C(ii,jj,:) + reshape(conv(reshape(A(ii,kk,:),[],1),reshape(B(kk,jj,:),[],1)),1,1,[]);
                end
            end
        end
    else
        error('Error: Check dimensions of polynomial matrices!');
    end
end

function A = parahermitian(R)
    % applies parahermitian that is the hermitian transpose + time-reversal
    A = flip(R,3);
    for ii=1:size(A,3)
        A(:,:,ii) = A(:,:,ii)';
    end
end



%% speech enhancement utility functions

function R = frame_stcov(x,T,W)
    % compute space-time covariance matrix accumulated over all frames
    [N, M] = size(x);
    R = zeros(M,M,2*T+1);
    framesize = T+1;
    for ff = 0:floor(N/framesize)-1
        xin = x(ff*framesize+1:(ff+1)*framesize,:);
        R = R + stcov(xin,W);
    end
    R = R./ff; % average across frames
end

function Rxx = stcov(xin,W)
    % computes the space-time covariance matrix using multi-channel frames
    [T,M] = size(xin);
    R = zeros(M,M,2*W+1);
    for ii=1:M
        for jj=1:M
            Rxx(ii,jj,:) = xcorr(xin(:,ii),xin(:,jj),W,'unbiased');
        end
    end
end

function  R = pm_recon(U,D)
    % reconstruct the polynomial matrix R = U^P * D * U using U and D
    R = pmdpdt(parahermitian(U),pmdpdt(D,U));
end
