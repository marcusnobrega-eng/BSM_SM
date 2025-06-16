function [D_Matrix] = Find_D_Matrix(f_Dir,coord_outlet,D_Matrix)
% Calculates the a matrix D_matrix where
% Each 1 represents an inflow into a cell and each -1, represents the outflow from the
% cell
% The (Inflow - Outflow) is given by D_matrix * q(t), where q(t)
% concatenates all flows into a row vector with dimention = size(F_dir)

dimensions = size(f_Dir);
row_out = coord_outlet(1); col_out = coord_outlet(2);
n_rows = dimensions(1);
n_cols = dimensions(2);
k = 0;
% z = f_Dir(:);
% D_Matrix(diag(z>0)) = -1;
for col = 1 : n_cols
    for row = 1 : n_rows
        k = k + 1;
        a = 0 ; b = 0 ; c = 0 ; d = 0 ; e = 0 ; f = 0 ; g = 0 ; h = 0 ; i = 0;
        if ((row == 1) && (col == 1))
            %a = e;
            %b = e;
            %c = e;
            %d = e;
            f = f_Dir(row, col + 1);
            %g = e;
            h = f_Dir(row + 1, col);
            i = f_Dir(row + 1, col + 1);
        elseif ((row == 1) && (col > 1) && (col < n_cols))
            %a = e;
            %b = e;
            %c = e;
            d = f_Dir(row, col - 1);
            f = f_Dir(row, col + 1);
            g = f_Dir(row + 1, col - 1);
            h = f_Dir(row + 1, col);
            i = f_Dir(row + 1, col + 1);
        elseif ((row == 1) && (col == n_cols))
            %a = e;
            %b = e;
            %c = e;
            d = f_Dir(row, col - 1);
            %f = e;
            g = f_Dir(row + 1, col - 1);
            h = f_Dir(row + 1, col);
            %i = e;
        elseif ((row > 1) && (row < n_rows) && (col == 1))
            %a = e;
            b = f_Dir(row - 1, col);
            c = f_Dir(row - 1, col + 1);
            %d = e;
            f = f_Dir(row, col + 1);
            %g = e;
            h = f_Dir(row + 1, col);
            i = f_Dir(row + 1, col + 1);
        elseif ((row > 1) && (row < n_rows) && (col == n_cols))
            a = f_Dir(row - 1, col - 1);
            b = f_Dir(row - 1, col);
            %c = e;
            d = f_Dir(row, col - 1);
            %f = e;
            g = f_Dir(row + 1, col - 1);
            h = f_Dir(row + 1, col);
            %i = e;
        elseif ((row == n_rows) && (col == 1))
            %a = e;
            b = f_Dir(row - 1, col);
            c = f_Dir(row - 1, col + 1);
            %d = e;
            f = f_Dir(row, col + 1);
            %g = e;
            %h = e;
            %i = e;
        elseif ((col > 1) && (col < n_cols) && (row == n_rows))
            a = f_Dir(row - 1, col - 1);
            b = f_Dir(row - 1, col);
            c = f_Dir(row - 1, col + 1);
            d = f_Dir(row, col - 1);
            f = f_Dir(row, col + 1);
            %g = e;
            %h = e;
            %i = e;
        elseif ((col == n_cols) && (row == n_rows))
            a = f_Dir(row - 1, col - 1);
            b = f_Dir(row - 1, col);
            %c = e;
            d = f_Dir(row, col - 1);
            %f = e;
            %g = e;
            %h = e;
            %i = e;
        else
            a = f_Dir(row - 1, col - 1);
            b = f_Dir(row - 1, col);
            c = f_Dir(row - 1, col + 1);
            d = f_Dir(row, col - 1);
            f = f_Dir(row, col + 1);
            g = f_Dir(row + 1, col - 1);
            h = f_Dir(row + 1, col);
            i = f_Dir(row + 1, col + 1);
        end

        if (a == 2) % Converting to 1
            D_Matrix(k,k - n_rows - 1) = 1;
        end
        if (b == 4)
            D_Matrix(k,k - 1) = 1;
        end
        if (c == 8)
            D_Matrix(k,k + n_rows - 1) = 1;
        end
        if (d == 1)
            D_Matrix(k,k - n_rows) = 1;
        end
        if (f == 16)
            D_Matrix(k,k + n_rows) = 1;
        end
        if (g == 128)
            D_Matrix(k,k - n_rows + 1) = 1;
        end
        if (h == 64)
            D_Matrix(k,k + 1) = 1;
        end
        if (i == 32)
            D_Matrix(k,k + n_rows + 1) = 1;
        end
%         if f_Dir(row,col) == 0
%             D_Matrix(k,k) = 0; % No outflow direction  
%         end
%         if isnan(f_Dir(row,col)) %%%% We are taking this out. It seems we
%         don't need it
%             D_Matrix(k,k) = 0; % Nan value
%         end
%         if f_Dir(row,col) > 0
%             D_Matrix(k,k) = -1; % Flow direction towards an dowstream cell
%         end        
        if row == row_out && col == col_out
            D_Matrix(k,k) = -1; % Outlet
        end
    end
end
end

