//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Duc Nguyen, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_1D
//##################################################################### 
#include <Grids/GRID_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void GRID_1D<T>::
Initialize(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid)
{
    m=m_input;xmin=xmin_input;xmax=xmax_input;
    if(MAC_grid){
        MAC_offset=(T).5;number_of_cells_x=max(m,0);
        if(m>0){dx=(xmax-xmin)/m;one_over_dx=1/dx;}
        else{dx=0;one_over_dx=0;}}
    else{ // regular grid
        MAC_offset=0;number_of_cells_x=max(m-1,0);
        if(number_of_cells_x){dx=(xmax-xmin)/(m-1);one_over_dx=1/dx;}
        else{dx=0;one_over_dx=0;}}
}
//#####################################################################
// Function Create_Grid_Given_Cell_Size
//#####################################################################
template<class T> GRID_1D<T> GRID_1D<T>::
Create_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create_Even_Sized_Grid_Given_Cell_Size
//#####################################################################
template<class T> GRID_1D<T> GRID_1D<T>::
Create_Even_Sized_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> template<class RW> void GRID_1D<T>::
Read(std::istream& input)
{
    Read_Binary<RW>(input,m);PHYSBAM_ASSERT(m>0);
    Read_Binary<RW>(input,xmin,xmax,MAC_offset);
    Initialize(m,xmin,xmax,MAC_offset!=0);
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class RW> void GRID_1D<T>::
Write(std::ostream& output) const
{
    Write_Binary<RW>(output,m,xmin,xmax,MAC_offset);
}
//#####################################################################
template class GRID_1D<float>;
template void GRID_1D<float>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_1D<float>::Write<float>(std::basic_ostream<char,std::char_traits<char> >&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_1D<double>;
template void GRID_1D<double>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_1D<double>::Write<double>(std::basic_ostream<char,std::char_traits<char> >&) const;
template void GRID_1D<double>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_1D<double>::Write<float>(std::basic_ostream<char,std::char_traits<char> >&) const;
template void GRID_1D<float>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_1D<float>::Write<double>(std::basic_ostream<char,std::char_traits<char> >&) const;
#endif
