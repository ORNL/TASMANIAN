function tsgWriteMatrix(filename, mat)
%
% tsgWriteMatrix(filename, mat)
%
% write a matrix to a text file
%
% the matrix [1 2 3 4; 5 6 7 8; 9 10 11 12;]
%
% is written as
%
% 3 4
% 1 2 3 4
% 5 6 7 8
% 9 10 11 12
%

filename = regexprep(filename, '\\ ', ' ');

if (prod(size(mat)) < 1000) % small matrix, use ascii format

    fid = fopen(filename, 'w');

    fprintf(fid, '%d  %d\n', size(mat, 1), size(mat, 2)); % load the number of points

    %format long;

    Ni = size(mat, 1);
    Nj = size(mat, 2);

    fmt = [''];

    for i = 1:Nj
        fmt = [fmt, '%2.20e '];
    end

    fmt = [fmt(1:end-1), '\n'];

    for i = 1:Ni
        fprintf(fid, fmt, mat(i,:));
    end

    fclose(fid);

else

    fid = fopen(filename, 'wb');

    Ni = size(mat, 1);
    Nj = size(mat, 2);

    fwrite(fid, ['TSG']);
    fwrite(fid, [Ni, Nj], 'integer*4');
    fwrite(fid, mat', 'double');

    %if (size(mat, 1) > 10)
    %    size(mat)
    %    mat(1:10, :)
    %end

    fclose(fid);

end

end
