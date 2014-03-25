function R = gDCA(filename, varargin)

% GDCA performs Gaussian Direct Coupling Analysis on a given FASTA file
%
% This function reads a filename in FASTA format and returns a ranking
% of residue pairs, as a matrix with 3 columns. Each column contains
% two indices i and j and a score. The rows are sorted with higer scores
% on top.
%
% The input file can be in plain text format or, if the provided
% MEX file was compiled, compressed via gzip.
%
% Example: GDCA('PF00014.fasta')
%
% Some optional parameters can be passed:
%
%   * pseudocount: the value of the pseudo-count parameter, between 0 and 1.
%                  the default is 0.8, which gives good results when the
%                  Frobenius norm score is used (see below); a good value for the
%                  Direct Information score is 0.2.
%   * theta: the value of the similarity threshold. By default it is 'auto',
%            which means it will be automatically computed (this takes additional
%            time); otherwise, a real value between 0 and 1 can be given.
%   * max_gap_fraction: maximum fraction of gap symbols in a sequence; sequences
%                       which exceed this threshold are discarded. The default
%                       value is 0.9.
%   * score: the scoring function to use. There are two possibilities, 'DI' for
%            the Direct Information, and 'frob' for the Frobenius norm. The
%            default is 'frob'.
%   * min_separation`: the minimum separation between residues in the output
%                      ranking. Must be >= 1. The default is 5.
%
% Example: GDCA('PF00014.fasta', 'pseudocount', 0.2, 'score', 'DI')
%
% This code accompanies the paper "Fast and accurate multivariate
% Gaussian modeling of protein families: Predicting residue contacts and
% protein-interaction partners" by Carlo Baldassi, Marco Zamparo, Christoph
% Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt and Andrea
% Pagnani, (2014) PLoS ONE 9(3): e92721. %doi:10.1371/journal.pone.0092721
%
% The code is released under the terms of the GNU General Public License
% version 3 or, at your option, any later version.

    p = inputParser;
    def_pseudocount = 0.8;
    def_max_gap_fraction = 0.9;
    def_score = 'frob';
    def_theta = 'auto';
    def_min_separation = 5;

    valid_scores = {'frob', 'DI'};
    valid_thetas = {'auto'};

    addRequired(p, 'filename', @ischar);
    addOptional(p, 'pseudocount', def_pseudocount, @(x) (isnumeric(x) && x >= 0 && x <=1));
    addOptional(p, 'max_gap_fraction', def_max_gap_fraction, @(x) (isnumeric(x) && x >= 0 && x <=1));
    addOptional(p, 'score', def_score, @(x) any(validatestring(x, valid_scores)));
    addOptional(p, 'theta', def_theta, ...
        @(x) ((isnumeric(x) && x >= 0 && x <= 1) || any(validatestring(x, valid_thetas))));
    addOptional(p, 'min_separation', def_min_separation, ...
        @(x) (isnumeric(x) && x >= 1));

    parse(p, filename, varargin{:});

    f = fopen(filename);
    if f == -1
        error('file not found: %s', filename);
    end
    fclose(f);

    pseudocount = p.Results.pseudocount;
    max_gap_fraction = p.Results.max_gap_fraction;
    score = validatestring(p.Results.score, valid_scores);
    theta = p.Results.theta;
    min_separation = p.Results.min_separation;

    path(path, 'modules');

    [N, ~, Z] = read_alignment_fasta(filename, max_gap_fraction);

    q = max(max(Z));

    [mJ, C] = compute_J(Z, pseudocount, theta);

    if strcmp(score, 'DI')
        S = compute_DI(mJ, C, N, q);
    else
        S = compute_frob(mJ, N, q);
    end

    S = correct_APC(S);

    R = compute_ranking(S, min_separation);
end

function [mJ, C] = compute_J(Z, pc, theta)

    [M, N] = size(Z);
    q = max(max(Z));

    [Pij_true, Pi_true, Wrew] = compute_true_frequencies(Z, M, N, q, theta);

    [Pij, Pi] = add_pseudocount(Pij_true, Pi_true, pc, N, q);

    C = compute_C(Pij, Pi);

    mJ = inv(C);

end

function [Pij, Pi] = add_pseudocount(Pij_true, Pi_true, pc, N, q)

    Pij = (1. - pc) * Pij_true + pc / q^2;
    Pi = (1. - pc) * Pi_true + pc / q;

    s = q - 1;

    i0 = 0;
    for i = 1 : N
	xr = i0 + (1:s);
	Pij(xr, xr) = (1. - pc) * Pij_true(xr, xr);
        for alpha = 1 : s
            x = i0 + alpha;
            Pij(x, x) = Pij(x, x) + pc / q;
        end
        i0 = i0 + s;
    end

end

function C = compute_C(Pij, Pi)
    C = Pij - Pi * Pi';
end

function iK = invsqrt(K)
    [U, L] = eig(K);
    iK = U * sqrt(L) * U';
end

function DI = compute_DI(mJ, Cmat, N, q)

    DI = zeros(N, N);
    s = q - 1;

    iKs = cell(N);

    rowi = 1;
    for i = 1 : N
        row = rowi : rowi + s - 1;
        rowi = rowi + s;

        Ki = Cmat(row, row);
        iKs{i} = invsqrt(Ki);
    end

    z = 0.5 * s * log(0.5);

    for i = 1 : N - 1
        row = (i-1)*s+1 : i*s;

        invsqrtKi = iKs{i};
        for j = i + 1 : N
            col = (j-1)*s+1 : j*s;

            invsqrtKj = iKs{j};
            mJij = mJ(row, col);
            if sum(sum(mJij)) ~= 0
                MM = invsqrtKi * mJij * invsqrtKj;
                V = MM * MM';
		[~, eigV] = eig(V);
		eigX = sqrt(1 + 4 * diag(eigV));
                DI(i,j) = z + 0.5 * sum(log(1 + eigX));
                DI(j,i) = DI(i,j);
            end
        end
    end
    for i = 1 : N
        DI(i,i) = 0;
    end

end

function F = compute_frob(mJ, N, q)

    F = zeros(N, N);
    s = q - 1;

    for i = 1 : N - 1
        row = (i-1)*s+1 : i*s;

        for j = i + 1 : N
            col = (j-1)*s+1 : j*s;

            mJij = mJ(row,col);
            mK = mJij - ones(s,1) * mean(mJij,1) - mean(mJij,2) * ones(1,s) + mean(mean(mJij));

            F(i,j) = norm(mK, 'fro');
            F(j,i) = F(i,j);
        end
    end

end

function S = correct_APC(S)
    N = size(S, 1);
    Si = sum(S, 1);
    Sj = sum(S, 2);
    Sa = sum(sum(S)) * (1 - 1/N);

    S = S - (Sj * Si) / Sa;
end

function R = compute_ranking(S, min_separation)
    N = size(S, 1);
    UR = zeros((N-min_separation)*(N-min_separation+1) / 2, 3);

    counter = 0;
    for i = 1 : N-min_separation
        for j = i + min_separation : N
            counter = counter + 1;
            UR(counter, 1) = i;
            UR(counter, 2) = j;
            UR(counter, 3) = S(i, j);
        end
    end

    [~,I] = sort(UR(:, 3), 'descend');
    R = UR(I,:);
end
