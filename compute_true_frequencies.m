function [Pij_true, Pi_true, Wrew] = compute_true_frequencies(Z, M, N, q, theta)

    try
        if isnumeric(theta)
            ntheta = theta;
        else
            ntheta = -1;
        end
        [Pij_true, Pi_true, Wrew] = compute_true_frequencies_mex(Z', M, N, q, ntheta);
    catch EM
        if isequal(EM.identifier, 'MATLAB:UndefinedFunction')
            warning('GaussDCA:fallback', 'compute_true_frequencies module not found, using slower fallback');
            [Pij_true, Pi_true, Wrew] = compute_true_frequencies_mat(Z, M, N, q, theta);
        else
            rethrow(EM);
        end
    end

end

function [Pij_true, Pi_true, W] = compute_true_frequencies_mat(Z, M, N, q, theta)

    if ~isnumeric(theta)
        theta = compute_theta(Z, M, N);
    end
    if theta > 0
        [W, Meff] = compute_reweighting(Z, M, N, theta);
    else
        W = ones(1, M);
        Meff = M;
    end

    [Pij_true, Pi_true] = compute_freqs(Z, W, q, Meff);

end

function theta = compute_theta(Z, M, N)
    theta = 0.38 * 0.32 / mean(1 - pdist(Z, 'hamm'));
    fprintf('theta=%f\n', theta);
end

function [W, Meff] = compute_reweighting(Z, M, N, theta)
    W = (1 ./ (1 + sum(squareform(pdist(Z, 'hamm') < (floor(N * theta)) / N))));

    Meff = sum(W);

    fprintf('M = %i N = %i Meff = %g\n', M, N, Meff);
end

function [Pij, Pi] = compute_freqs(Z, W, q, Meff)
    [M, N] = size(Z);
    s = q - 1;

    Ns = N * s;

    Pij = zeros(Ns, Ns);
    Pi = zeros(Ns, 1);

    i0 = 0;
    for i = 1:N
        for k = 1:M
            a = Z(k,i);
            if a ~= q
                Pi(i0 + a) = Pi(i0 + a) + W(k);
            end
        end
        i0 = i0 + s;
    end
    Pi = Pi ./ Meff;

    i0 = 0;
    for i = 1:N
        j0 = i0;
        for j = i:N
            for k = 1:M
                a = Z(k,i);
                b = Z(k,j);
                if a ~= q && b ~= q
                    Pij(i0+a, j0+b) = Pij(i0+a, j0+b) + W(k);
                end
            end
            j0 = j0 + s;
        end
        i0 = i0 + s;
    end
    for i = 1:Ns
        Pij(i,i) = Pij(i,i) ./ Meff;
        for j = i+1:Ns
            Pij(i,j) = Pij(i,j) ./ Meff;
            Pij(j,i) = Pij(i,j);
        end
    end
end
