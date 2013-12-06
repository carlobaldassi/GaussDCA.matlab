function [N, M, Z] = read_alignment_fasta(filename, max_gap_fraction)
    try
        [N, M, Z] = read_alignment_fasta_mex(filename, max_gap_fraction);
    catch EM
        if isequal(EM.identifier, 'MATLAB:UndefinedFunction')
            warning('GaussDCA:fallback', 'read_alignment_fasta module not found, using slower fallback');
            [N, M, Z] = read_alignment_fasta_mat(filename, max_gap_fraction);
        else
            rethrow(EM);
        end
    end
    Z = Z(1:M,:);
end

function [N, M, Z] = read_alignment_fasta_mat(filename, max_gap_fraction)
    align_full = fastaread(filename);
    M = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Z = zeros(M,N);
    zcounter = 0;
    for i = 1:M
        counter = 0;
        zs = zeros(1,N);
        numgaps = 0;
        for j = 1:length(ind)
            if ind(j)
                counter = counter + 1;
                l = align_full(i).Sequence(j);
                numgaps = numgaps + isgap(l);
                zs(counter) = letter2number(l);
            end
        end
        if numgaps / N > max_gap_fraction
            continue;
        end
        zcounter = zcounter + 1;
        Z(zcounter,:) = zs;
    end
    M = zcounter;
end

function x = letter2number(a)
    switch a
        case 'A'
            x = 1;
        case 'C'
            x = 2;
        case 'D'
            x = 3;
        case 'E'
            x = 4;
        case 'F'
            x = 5;
        case 'G'
            x = 6;
        case 'H'
            x = 7;
        case 'I'
            x = 8;
        case 'K'
            x = 9;
        case 'L'
            x = 10;
        case 'M'
            x = 11;
        case 'N'
            x = 12;
        case 'P'
            x = 13;
        case 'Q'
            x = 14;
        case 'R'
            x = 15;
        case 'S'
            x = 16;
        case 'T'
            x = 17;
        case 'V'
            x = 18;
        case 'W'
            x = 19;
        case 'Y'
            x = 20;
        case '-'
            x = 21;
        otherwise
            x = 21;
    end
end

function tf = isgap(l)
    tf = (l == '-');
end
