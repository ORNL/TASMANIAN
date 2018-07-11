function tsgWriteMatrix(filename, mat)
%
% tsgWriteMatrix(filename, mat)
%
% write a matrix to a text file
% 
% Both the real matrix [1 2 3 4; 5 6 7 8; 9 10 11 12] and the complex
% matrix [1+2i 3+4i; 5+6i 7+8i; 9+10i 11+12i]
%
% are written as
%
% 3 4
% 1 2 3 4
% 5 6 7 8
% 9 10 11 12
%

bIsComplex = ~isreal(mat);

if (prod(size(mat)) < 1000) % small matrix, use ascii format

    fid = fopen(filename, 'w');

    Ni = size(mat, 1);

    if(bIsComplex)
        Nj = 2 * size(mat, 2);
    else
        Nj = size(mat, 2);
    end

    fprintf(fid, '%d  %d\n', Ni, Nj); % load the number of points

    %format long;

    fmt = [''];

    for i = 1:Nj
        if (~bIsComplex)
            fmt = [fmt, '%2.20e '];
        else
            fmt = [fmt, '%2.20e %2.20e '];
        end
    end

    fmt = [fmt(1:end-1), '\n'];

    if (~bIsComplex)
        for i = 1:Ni
            fprintf(fid, fmt, mat(i,:));
        end
    else
        pmat = zeros(Ni, Nj);
        pmat(:,1:2:Nj) = real(mat);
        pmat(:,2:2:Nj) = imag(mat);
        for i = 1:Ni
            fprintf(fid, fmt, pmat(i,:));
        end
    end

    fclose(fid);

else

    fid = fopen(filename, 'wb');
    
    Ni = size(mat, 1);

    if (bIsComplex)
        Nj = 2 * size(mat, 2);
    else
        Nj = size(mat, 2);
    end
    
    fwrite(fid, ['TSG']);
    fwrite(fid, [Ni, Nj], 'integer*4');
    if (~bIsComplex)
        fwrite(fid, mat', 'double');
    else
        pmat = zeros(Ni,Nj);
        pmat(:,1:2:Nj) = real(mat);
        pmat(:,2:2:Nj) = imag(mat);
        fwrite(fid, pmat', 'double');
    end
    
    %if (size(mat, 1) > 10)
    %    size(mat)
    %    mat(1:10, :)
    %end
    
    fclose(fid);

end

end
