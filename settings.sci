function [y] = settings(a)

//    a = gca();
    a.font_size = 2;         // Increase tick number size
    a.x_label.font_size = 3; // Increase X label size
    a.y_label.font_size = 3; // Increase Y label size
    xgrid;
    y = 0;

endfunction
