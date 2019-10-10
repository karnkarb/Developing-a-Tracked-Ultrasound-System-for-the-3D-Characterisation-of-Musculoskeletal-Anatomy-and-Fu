function [sorted_cordinates] = sorting_lines(sort_line,coordinates_3)

    if sort_line(1) == coordinates_3(1,2) && sort_line(2) == coordinates_3(4,2)
                        X3 = coordinates_3(1, :);
                        X2 = coordinates_3(2, :);
                        X1 = coordinates_3(3, :);
                        X6 = coordinates_3(4, :);
                        X5 = coordinates_3(5, :);
                        X4 = coordinates_3(6, :);
                        X9 = coordinates_3(7, :);
                        X8 = coordinates_3(8, :);
                        X7 = coordinates_3(9, :);
    elseif sort_line(1) == coordinates_3(1,2) && sort_line(2) == coordinates_3(7,2)
                        X3 = coordinates_3(1, :);
                        X2 = coordinates_3(2, :);
                        X1 = coordinates_3(3, :);
                        X9 = coordinates_3(4, :);
                        X8 = coordinates_3(5, :);
                        X7 = coordinates_3(6, :);
                        X6 = coordinates_3(7, :);
                        X5 = coordinates_3(8, :);
                        X4 = coordinates_3(9, :);
    elseif sort_line(1) == coordinates_3(4,2) && sort_line(2) == coordinates_3(1,2)
                        X6 = coordinates_3(1, :);
                        X5 = coordinates_3(2, :);
                        X4 = coordinates_3(3, :);
                        X3 = coordinates_3(4, :);
                        X2 = coordinates_3(5, :);
                        X1 = coordinates_3(6, :);
                        X9 = coordinates_3(7, :);
                        X8 = coordinates_3(8, :);
                        X7 = coordinates_3(9, :);    
    elseif sort_line(1) == coordinates_3(4,2) && sort_line(2) == coordinates_3(7,2)
                        X9 = coordinates_3(1, :);
                        X8 = coordinates_3(2, :);
                        X7 = coordinates_3(3, :);
                        X3 = coordinates_3(4, :);
                        X2 = coordinates_3(5, :);
                        X1 = coordinates_3(6, :);
                        X6 = coordinates_3(7, :);
                        X5 = coordinates_3(8, :);
                        X4 = coordinates_3(9, :);   
    elseif sort_line(3) == coordinates_3(7,2) && sort_line(2) == coordinates_3(1,2)
                        X6 = coordinates_3(1, :);
                        X5 = coordinates_3(2, :);
                        X4 = coordinates_3(3, :);
                        X9 = coordinates_3(4, :);
                        X8 = coordinates_3(5, :);
                        X7 = coordinates_3(6, :);
                        X3 = coordinates_3(7, :);
                        X2 = coordinates_3(8, :);
                        X1 = coordinates_3(9, :);                    
    else
                        X6 = coordinates_3(1, :);
                        X5 = coordinates_3(2, :);
                        X4 = coordinates_3(3, :);
                        X9 = coordinates_3(4, :);
                        X8 = coordinates_3(5, :);
                        X7 = coordinates_3(6, :);
                        X3 = coordinates_3(7, :);
                        X2 = coordinates_3(8, :);
                        X1 = coordinates_3(9, :);
    end
    sorted_cordinates = [X1;X2;X3;X4;X5;X6;X7;X8;X9];
end