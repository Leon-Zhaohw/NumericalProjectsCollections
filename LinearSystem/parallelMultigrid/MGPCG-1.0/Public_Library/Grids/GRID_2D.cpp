//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_2D
//##################################################################### 
#include <Grids/GRID_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void GRID_2D<T>::
Initialize(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid)
{
    m=m_input;n=n_input;xmin=xmin_input;xmax=xmax_input;ymin=ymin_input;ymax=ymax_input;
    if(MAC_grid){
        MAC_offset=(T).5;number_of_cells_x=max(m,0);number_of_cells_y=max(n,0);
        if(m>0&&n>0){
            dx=(xmax-xmin)/m;dy=(ymax-ymin)/n;one_over_dx=1/dx;one_over_dy=1/dy;min_dx_dy=min(dx,dy);max_dx_dy=max(dx,dy);diagonal_length=sqrt(sqr(dx)+sqr(dy));}
        else{dx=0;dy=0;one_over_dx=0;one_over_dy=0;min_dx_dy=0;max_dx_dy=0;diagonal_length=0;}}
    else{ // regular grid
        MAC_offset=0;number_of_cells_x=max(m-1,0);number_of_cells_y=max(n-1,0);
        if(m>1&&n>1){
            dx=(xmax-xmin)/(m-1);dy=(ymax-ymin)/(n-1);one_over_dx=1/dx;one_over_dy=1/dy;min_dx_dy=min(dx,dy);max_dx_dy=max(dx,dy);diagonal_length=sqrt(sqr(dx)+sqr(dy));}
        else{dx=0;dy=0;one_over_dx=0;one_over_dy=0;min_dx_dy=0;max_dx_dy=0;diagonal_length=0;}}
}
//#####################################################################
// Function Create_Grid_Given_Cell_Size
//#####################################################################
template<class T> GRID_2D<T> GRID_2D<T>::
Create_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create_Even_Sized_Grid_Given_Cell_Size
//#####################################################################
template<class T> GRID_2D<T> GRID_2D<T>::
Create_Even_Sized_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> template<class RW> void GRID_2D<T>::
Read(std::istream& input)
{
    Read_Binary<RW>(input,m,n);PHYSBAM_ASSERT(m>0 && n>0);
    Read_Binary<RW>(input,xmin,xmax,ymin,ymax,MAC_offset);
    Initialize(m,n,xmin,xmax,ymin,ymax,MAC_offset!=0);
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class RW> void GRID_2D<T>::
Write(std::ostream& output) const
{
    Write_Binary<RW>(output,m,n,xmin,xmax,ymin,ymax,MAC_offset);
}
//#####################################################################
template class GRID_2D<float>;
template void GRID_2D<float>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_2D<float>::Write<float>(std::basic_ostream<char,std::char_traits<char> >&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_2D<double>;
template void GRID_2D<double>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_2D<double>::Write<double>(std::basic_ostream<char,std::char_traits<char> >&) const;
template void GRID_2D<double>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_2D<double>::Write<float>(std::basic_ostream<char,std::char_traits<char> >&) const;
template void GRID_2D<float>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_2D<float>::Write<double>(std::basic_ostream<char,std::char_traits<char> >&) const;
#endif
