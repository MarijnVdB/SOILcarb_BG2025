function [YY,Yu] = calculate_YY_Yu(colNum, NU, NC, n, X_in, M)

    % YY and Yu are reconstructed
    % Yu are the results for the unconditional run, in an array
    Yu = X_in(1:NU,colNum);
    
    % YY are the results of the conditional runs, in a cell array
    YY = cell(M,n);
    % YY is filled
    counter = 1; % To keep track of the rows where the values have to be stored
    for i = 1:n % Columns
        for j = 1:M % Rows
            startRow = NU + NC*(counter-1) + 1;
            endRow = NU + (counter*NC);
            YY{j,i} = X_in(startRow:endRow, colNum);
            
%             YY{j,i} = X_in((counter*NU)+1:(NU*(counter+1)), colNum);
            counter = counter + 1;
        end
    end

end

