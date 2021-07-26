//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_3D
//##################################################################### 
#include <Grids/GRID_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void GRID_3D<T>::
Initialize(const int m_input,const int n_input,const int mn_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const T zmin_input,const T zmax_input,
    const bool MAC_grid)
{
    m=m_input;n=n_input;mn=mn_input;xmin=xmin_input;xmax=xmax_input;ymin=ymin_input;ymax=ymax_input;zmin=zmin_input;zmax=zmax_input;
    if(MAC_grid){
        MAC_offset=(T).5;number_of_cells_x=max(m,0);number_of_cells_y=max(n,0);number_of_cells_z=max(mn,0);
        if(m>0&&n>0&&mn>0){
            dx=(xmax-xmin)/m;dy=(ymax-ymin)/n;dz=(zmax-zmin)/mn;one_over_dx=1/dx;one_over_dy=1/dy;one_over_dz=1/dz;min_dx_dy_dz=min(dx,dy,dz);max_dx_dy_dz=max(dx,dy,dz);
            diagonal_length=sqrt(sqr(dx)+sqr(dy)+sqr(dz));}
        else{dx=0;dy=0;dz=0;one_over_dx=0;one_over_dy=0;one_over_dz=0;min_dx_dy_dz=0;max_dx_dy_dz=0;diagonal_length=0;}}
    else{ // regular grid
        MAC_offset=0;number_of_cells_x=max(m-1,0);number_of_cells_y=max(n-1,0);number_of_cells_z=max(mn-1,0);
        if(m>1&&n>1&&mn>1){
            dx=(xmax-xmin)/(m-1);dy=(ymax-ymin)/(n-1);dz=(zmax-zmin)/(mn-1);one_over_dx=1/dx;one_over_dy=1/dy;one_over_dz=1/dz;min_dx_dy_dz=min(dx,dy,dz);max_dx_dy_dz=max(dx,dy,dz);
            diagonal_length=sqrt(sqr(dx)+sqr(dy)+sqr(dz));}
        else{dx=0;dy=0;dz=0;one_over_dx=0;one_over_dy=0;one_over_dz=0;min_dx_dy_dz=0;max_dx_dy_dz=0;diagonal_length=0;}}
}
//#####################################################################
// Function Create_Grid_Given_Cell_Size
//#####################################################################
template<class T> GRID_3D<T> GRID_3D<T>::
Create_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Create_Even_Sized_Grid_Given_Cell_Size
//#####################################################################
template<class T> GRID_3D<T> GRID_3D<T>::
Create_Even_Sized_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> template<class RW> void GRID_3D<T>::
Read(std::istream& input)
{
    Read_Binary<RW>(input,m,n,mn);PHYSBAM_ASSERT(m>0 && n>0 && mn>0);
    Read_Binary<RW>(input,xmin,xmax,ymin,ymax,zmin,zmax,MAC_offset);
    Initialize(m,n,mn,xmin,xmax,ymin,ymax,zmin,zmax,MAC_offset!=0);
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> template<class RW> void GRID_3D<T>::
Write(std::ostream& output) const
{
    Write_Binary<RW>(output,m,n,mn);Write_Binary<RW>(output,xmin,xmax,ymin,ymax,zmin,zmax);Write_Binary<RW>(output,MAC_offset);
}
//#####################################################################
template class GRID_3D<float>;
template void GRID_3D<float>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_3D<float>::Write<float>(std::basic_ostream<char,std::char_traits<char> >&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_3D<double>;
template void GRID_3D<double>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_3D<double>::Write<double>(std::basic_ostream<char,std::char_traits<char> >&) const;
template void GRID_3D<double>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_3D<double>::Write<float>(std::basic_ostream<char,std::char_traits<char> >&) const;
template void GRID_3D<float>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
template void GRID_3D<float>::Write<double>(std::basic_ostream<char,std::char_traits<char> >&) const;
#endif
