function [dim_out,width_out,height_out]=dfig(n_row,m_column,text_width)
    
% dimension of the figure

fdim=(text_width/m_column);
dim_out=fdim;
% dimensions for the width and height parameters in  multiplot
width_out=(text_width/m_column);
height_out=dim_out*n_row;

end