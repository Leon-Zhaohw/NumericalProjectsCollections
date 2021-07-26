//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_2D
//#####################################################################
// A number of functions (e.g. the Clamp functions) assume the grid indexing starts at 1.  This way we can use truncation rather than floor because floor is really slow.
//#####################################################################
#ifndef __GRID_2D__
#define __GRID_2D__

#include <Grids/GRID_1D.h>
#include <Grids/GRID_POLICY.h>
#include <Arrays/ARRAYS_FORWARD.h>
#include <Geometry/BOX.h>
namespace PhysBAM{

template<class T_GRID> class LEVELSET_2D;
template<class T_GRID> class INCOMPRESSIBLE_2D;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_2D_HELPER;
class FLOOD_FILL_2D;

template<class T>
class GRID_2D:public POLICY_UNIFORM<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
public:
    enum REGION {WHOLE_REGION,GHOST_REGION,BOUNDARY_REGION,INTERIOR_REGION,BOUNDARY_INTERIOR_REGION}; // for iterators
    static const int dimension=2;
    static const int number_of_cells_per_block=4;
    static const int number_of_faces_per_block=12;
    static const int number_of_incident_faces_per_block=2;
    static const int number_of_nodes_per_face=2;
    static const int number_of_nodes_per_cell=4;
    static const int number_of_cells_per_node=4;
    static const int number_of_neighbors_per_node=4;
    static const int number_of_neighbors_per_cell=4;
    static const int number_of_faces_per_cell=4;
    static const int number_of_one_ring_neighbors_per_cell=8;
    typedef GRID_2D UNIFORM_GRID;
    typedef UNIFORM_GRID_ITERATOR_NODE_2D<T> NODE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_CELL_2D<T> CELL_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_FACE_2D<T> FACE_ITERATOR;
    typedef LEVELSET_2D<GRID_2D> LEVELSET;
    typedef INCOMPRESSIBLE_2D<GRID_2D> INCOMPRESSIBLE;
    typedef FLOOD_FILL_2D FLOOD_FILL;
    typedef LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID_2D> LINEAR_INTERPOLATION_MAC_HELPER;

    int m,n; // # of points: x and y direction
    T xmin,xmax,ymin,ymax; // left and right wall, bottom and top wall
    T dx,dy; // cell sizes
    T one_over_dx,one_over_dy,min_dx_dy,max_dx_dy,diagonal_length;
    T MAC_offset; // 0 for a regular grid and .5 for a MAC grid
    int number_of_cells_x,number_of_cells_y; // not saved to file

    GRID_2D()
    {
        Initialize(TV_INT(),BOX<TV>::Unit_Box());
    }

    GRID_2D(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid=false)
    {
        Initialize(m_input,n_input,xmin_input,xmax_input,ymin_input,ymax_input,MAC_grid);
    }

    GRID_2D(const int m_input,const int n_input,const BOX<TV>& box,const bool MAC_grid=false)
    {
        Initialize(m_input,n_input,box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,MAC_grid);
    }

    GRID_2D(const TV_INT& counts,const BOX<TV>& box,const bool MAC_grid=false)
    {
        Initialize(counts,box,MAC_grid);
    }

    GRID_2D(const GRID_2D<T>& grid_input)
    {
        Initialize(grid_input.Counts(),grid_input.Domain(),grid_input.Is_MAC_Grid());
    }

    template<class T2>
    GRID_2D(const GRID_2D<T2>& grid_input)
    {
        Initialize(grid_input.Counts(),BOX<TV>(grid_input.Domain()),grid_input.Is_MAC_Grid());
    }

    GRID_2D<T>& operator=(const GRID_2D<T>& grid_input)
    {if(this!=&grid_input) Initialize(grid_input.Counts(),grid_input.Domain(),grid_input.Is_MAC_Grid());return *this;}

    void Set_To_Double_Resolution_Grid(const GRID_2D<T>& grid)
    {Initialize(2*grid.Numbers_Of_Cells()+!grid.Is_MAC_Grid(),grid.Domain(),grid.Is_MAC_Grid());}

    void Initialize(const TV_INT& counts,const BOX<TV>& box,const bool MAC_grid=false)
    {Initialize(counts.x,counts.y,box,MAC_grid);}

    void Initialize(const int m_input,const int n_input,const BOX<TV>& box,const bool MAC_grid=false)
    {Initialize(m_input,n_input,box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,MAC_grid);}

    bool Is_MAC_Grid() const
    {return MAC_offset==.5;}

    bool Is_Isotropic(const T tolerance=1e-4) const
    {return abs(dx-dy)<=tolerance*max_dx_dy;}

    TV DX() const
    {return TV(dx,dy);}

    TV One_Over_DX() const
    {return TV(one_over_dx,one_over_dy);}

    TV_INT Numbers_Of_Cells() const
    {return TV_INT(number_of_cells_x,number_of_cells_y);}

    TV_INT Numbers_Of_Nodes() const
    {return Numbers_Of_Cells()+1;}

    T Minimum_Edge_Length() const
    {return min_dx_dy;}

    T Maximum_Edge_Length() const
    {return max_dx_dy;}

    T Cell_Diagonal() const
    {return diagonal_length;}

    T Cell_Size() const
    {return dx*dy;}

    T Face_Size(const int axis) const
    {assert(1<=axis&&axis<=2);return DX()[3-axis];}

    TV Xmin() const
    {return TV(xmin,ymin);}

    TV Xmax() const
    {return TV(xmax,ymax);}

    T x(const int i) const // x at the i=1 to i=m grid points
    {return xmin+(i-(T)1+MAC_offset)*dx;}

    T y(const int j) const // y at the j=1 to j=n grid points
    {return ymin+(j-(T)1+MAC_offset)*dy;}

    T x_plus_half(const int i) const
    {return xmin+(i-(T).5+MAC_offset)*dx;}

    T x_minus_half(const int i) const
    {return xmin+(i-(T)1.5+MAC_offset)*dx;}

    T y_plus_half(const int j) const
    {return ymin+(j-(T).5+MAC_offset)*dy;}

    T y_minus_half(const int j) const
    {return ymin+(j-(T)1.5+MAC_offset)*dy;}

    TV X(const int i,const int j) const
    {return TV(x(i),y(j));}

    TV X(const TV_INT& index) const
    {return TV(x(index.x),y(index.y));}

    TV Node(const int i,const int j) const
    {return TV(xmin+(i-(T)1)*dx,ymin+(j-(T)1)*dy);}

    TV Node(const TV_INT& index) const
    {return TV(xmin+(index.x-(T)1)*dx,ymin+(index.y-(T)1)*dy);}

    TV Center(const int i,const int j) const
    {return TV(xmin+(i-(T).5)*dx,ymin+(j-(T).5)*dy);}

    TV Center(const TV_INT& index) const
    {return TV(xmin+(index.x-(T).5)*dx,ymin+(index.y-(T).5)*dy);}

    TV X_Face(const int i,const int j) const
    {return TV(xmin+(i-(T)1)*dx,ymin+(j-(T).5)*dy);}

    TV X_Face(const TV_INT& index) const
    {return TV(xmin+(index.x-(T)1)*dx,ymin+(index.y-(T).5)*dy);}

    TV Y_Face(const int i,const int j) const
    {return TV(xmin+(i-(T).5)*dx,ymin+(j-(T)1)*dy);}

    TV Y_Face(const TV_INT& index) const
    {return TV(xmin+(index.x-(T).5)*dx,ymin+(index.y-(T)1)*dy);}

    TV Face(const int axis,const TV_INT& index) const
    {switch(axis){
      case 1:return X_Face(index);
      default:assert(axis==2);return Y_Face(index);}}

    BOX<TV> Face_Domain(const int axis,const TV_INT& index,const T thickness_over_two=0) const
    {TV dimensions=(T).5*DX();dimensions[axis]=thickness_over_two;BOX<TV> domain(Face(axis,index));domain.Change_Size(dimensions);return domain;}

    static TV_INT Node_Cell_Index(const TV_INT& node,const int cell)
    {static const TV_INT cell_from_node_offset[4]={TV_INT(-1,-1),TV_INT(0,-1),TV_INT(-1,0),TV_INT(0,0)};
    assert(1<=cell&&cell<=4);return node+cell_from_node_offset[cell-1];}

    static TV_INT Face_Node_Index(const int axis,const TV_INT& face_index,const int node)
    {assert(1<=node&&node<=2);TV_INT index=face_index;index[3-axis]+=node-1;return index;}

    static TV_INT Node_Face_Index(const int axis,const TV_INT& node_index,const int face)
    {assert(1<=face&&face<=2);TV_INT index=node_index;index[3-axis]+=face-2;return index;}

    void Face_Corner_To_Opposite_Corner_Vectors(const int axis,TV vectors[2])
    {vectors[1]=DX();vectors[1][axis]=0;vectors[0]=(T)-1*vectors[1];}

    T Face_Corner_To_Opposite_Corner_Length(const int axis) const
    {return DX()[3-axis];}

    void Index(const TV& location,int& i,int& j) const // returns the left and bottom indices
    {i=(int)floor((location.x-xmin)*one_over_dx+1-MAC_offset); // note that "floor" is expensive
    j=(int)floor((location.y-ymin)*one_over_dy+1-MAC_offset);}

    void Cell(const TV& location,int& i,int& j,const int number_of_ghost_cells) const // returns the left and bottom
    {int number_of_ghost_cells_plus_one=number_of_ghost_cells+1;
    i=(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells_plus_one)-number_of_ghost_cells;
    j=(int)((location.y-ymin)*one_over_dy+number_of_ghost_cells_plus_one)-number_of_ghost_cells;}

    TV_INT Cell(const TV& location,const int number_of_ghost_cells) const // returns the left and bottom
    {TV_INT index;Cell(location,index.x,index.y,number_of_ghost_cells);return index;}

    BOX<TV> Cell_Domain(const int i,const int j) const
    {return BOX<TV>(xmin+(i-(T)1)*dx,xmin+i*dx,ymin+(j-(T)1)*dy,ymin+j*dy);}

    BOX<TV> Cell_Domain(const TV_INT& index) const
    {return BOX<TV>(xmin+(index.x-(T)1)*dx,xmin+index.x*dx,ymin+(index.y-(T)1)*dy,ymin+index.y*dy);}

    static void Cells_Touching_Face(const int axis,const TV_INT& face_index,TV_INT& cell1,TV_INT& cell2)
    {cell2=face_index;cell1=face_index;cell1[axis]-=1;}

    TV Clamp(const TV& location) const
    {return TV(clamp(location.x,xmin,xmax),clamp(location.y,ymin,ymax));}

    TV Clamp_Component(const TV& location,const int component) const
    {assert(1<=component&&component<=2);
    if(component==1)return TV(clamp(location.x,xmin,xmax),location.y);
    else return TV(location.x,clamp(location.y,ymin,ymax));}

    TV Clamp(const TV& location,int number_of_ghost_cells) const // clamps to the grid (with ghost cells)
    {T extra_x=number_of_ghost_cells*dx,extra_y=number_of_ghost_cells*dy;
    return TV(clamp(location.x,xmin-extra_x,xmax+extra_x),clamp(location.y,ymin-extra_y,ymax+extra_y));}

    TV_INT Clamp_Min(const TV_INT& index) const
    {return TV_INT(clamp_min(index.x,1),clamp_min(index.y,1));}

    TV_INT Clamp_Max(const TV_INT& index) const
    {return TV_INT(clamp_max(index.x,m),clamp_max(index.y,n));}

    void Clamp(int& i,int& j) const
    {i=clamp(i,1,m);j=clamp(j,1,n);}

    void Clamp(TV_INT& index) const
    {index.x=clamp(index.x,1,m);index.y=clamp(index.y,1,n);}

    void Clamp_End_Minus_One(int& i,int& j) const
    {i=clamp(i,1,m-1);j=clamp(j,1,n-1);}

    void Clamped_Index(const TV& location,int& i,int& j) const
    {i=min(m,1+max(0,(int)((location.x-xmin)*one_over_dx-MAC_offset)));
    j=min(n,1+max(0,(int)((location.y-ymin)*one_over_dy-MAC_offset)));}

    TV_INT Clamped_Index(const TV& location) const
    {TV_INT index;Clamped_Index(location,index.x,index.y);return index;}

    void Clamped_Index_End_Minus_One(const TV& location,int& i,int& j) const
    {i=min(m-1,1+max(0,(int)((location.x-xmin)*one_over_dx-MAC_offset)));
    j=min(n-1,1+max(0,(int)((location.y-ymin)*one_over_dy-MAC_offset)));}

    TV_INT Clamped_Index_End_Minus_One(const TV& location) const
    {TV_INT index;Clamped_Index_End_Minus_One(location,index.x,index.y);return index;}

    void Clamped_Index_End_Minus_One(const TV& location,int& i,int& j,const int number_of_ghost_cells) const
    {i=min(m+number_of_ghost_cells-1,1-number_of_ghost_cells+max(0,(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells-MAC_offset)));
    j=min(n+number_of_ghost_cells-1,1-number_of_ghost_cells+max(0,(int)((location.y-ymin)*one_over_dy+number_of_ghost_cells-MAC_offset)));}

    TV_INT Clamped_Index_End_Minus_One(const TV& location,const int number_of_ghost_cells) const
    {TV_INT index;Clamped_Index_End_Minus_One(location,index.x,index.y,number_of_ghost_cells);return index;}

    void Clamp_To_Cell(const TV& location,int& i,int& j) const
    {i=min(number_of_cells_x,1+max(0,(int)((location.x-xmin)*one_over_dx)));
    j=min(number_of_cells_y,1+max(0,(int)((location.y-ymin)*one_over_dy)));}

    TV_INT Clamp_To_Cell(const TV& location) const
    {TV_INT index;Clamp_To_Cell(location,index.x,index.y);return index;}

    void Clamp_To_Cell(const TV& location,int& i,int& j,const int number_of_ghost_cells) const
    {i=min(number_of_cells_x+number_of_ghost_cells,1-number_of_ghost_cells+max(0,(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells)));
    j=min(number_of_cells_y+number_of_ghost_cells,1-number_of_ghost_cells+max(0,(int)((location.y-ymin)*one_over_dy+number_of_ghost_cells)));}

    TV_INT Clamp_To_Cell(const TV& location,const int number_of_ghost_cells) const
    {TV_INT index;Clamp_To_Cell(location,index.x,index.y,number_of_ghost_cells);return index;}

    BOX<TV_INT> Clamp_To_Cell(const BOX<TV>& box,const int number_of_ghost_cells) const
    {return BOX<TV_INT>(Clamp_To_Cell(box.Minimum_Corner(),number_of_ghost_cells),Clamp_To_Cell(box.Maximum_Corner(),number_of_ghost_cells));}

    TV_INT Block_Index(const TV& X,const int number_of_ghost_cells) const // index of node at center of block
    {assert(Is_MAC_Grid());return Clamped_Index_End_Minus_One(X,number_of_ghost_cells)+TV_INT::All_Ones_Vector();}

    void Closest_Node(const TV& location,int& i,int& j) const
    {i=min(m,1+max(0,(int)((location.x-xmin)*one_over_dx+(T).5)));
    j=min(n,1+max(0,(int)((location.y-ymin)*one_over_dy+(T).5)));}

    TV_INT Closest_Node(const TV& location) const
    {TV_INT index;Closest_Node(location,index.x,index.y);return index;}

    BOX<TV> Domain() const
    {return BOX<TV>(xmin,xmax,ymin,ymax);}

    BOX<TV> Ghost_Domain(const int number_of_ghost_cells) const
    {return BOX<TV>(xmin-dx*number_of_ghost_cells,xmax+dx*number_of_ghost_cells,ymin-dy*number_of_ghost_cells,ymax+dy*number_of_ghost_cells);}

    TV_INT Counts() const
    {return TV_INT(m,n);}

    BOX<TV_INT> Domain_Indices() const
    {return BOX<TV_INT>(1,m,1,n);}

    BOX<TV_INT> Node_Indices(const int ghost_cells=0) const
    {return BOX<TV_INT>(TV_INT()+1,Numbers_Of_Nodes()).Thickened(ghost_cells);}

    BOX<TV_INT> Cell_Indices(const int ghost_cells=0) const
    {return BOX<TV_INT>(TV_INT()+1,Numbers_Of_Cells()).Thickened(ghost_cells);}

    BOX<TV_INT> Block_Indices(const int ghost_cells=0) const
    {return Node_Indices(ghost_cells);}

    VECTOR<BOX<TV_INT>,2> Face_Indices(const int ghost_cells=0) const
    {return VECTOR<BOX<TV_INT>,2>(Get_X_Face_Grid().Node_Indices(ghost_cells),Get_Y_Face_Grid().Node_Indices(ghost_cells));}

    bool Outside(const TV& location) const
    {return location.x < xmin || location.x > xmax || location.y < ymin || location.y > ymax;}

    GRID_2D<T> Get_MAC_Grid() const
    {return GRID_2D<T>(Numbers_Of_Cells(),Domain(),true);}

    GRID_2D<T> Get_Regular_Grid() const
    {return GRID_2D<T>(Numbers_Of_Nodes(),Domain(),false);}

    GRID_2D<T> Get_X_Face_Grid() const
    {return GRID_2D<T>(number_of_cells_x+1,max(0,number_of_cells_y),xmin,xmax,ymin+(T).5*dy,ymax-(T).5*dy);}

    GRID_2D<T> Get_Y_Face_Grid() const
    {return GRID_2D<T>(max(0,number_of_cells_x),number_of_cells_y+1,xmin+(T).5*dx,xmax-(T).5*dx,ymin,ymax);}

    GRID_2D<T> Get_Face_Grid(const int axis) const
    {switch(axis){case 1:return Get_X_Face_Grid();default:assert(axis==2);return Get_Y_Face_Grid();}}

    GRID_2D<T> Get_Regular_Grid_At_MAC_Positions() const
    {assert(Is_MAC_Grid());
    return GRID_2D<T>(m,n,xmin+(T).5*dx,xmax-(T).5*dx,ymin+(T).5*dy,ymax-(T).5*dy);}

    GRID_2D<T> Get_MAC_Grid_At_Regular_Positions() const
    {assert(!Is_MAC_Grid());return GRID_2D<T>(m,n,xmin-(T).5*dx,xmax+(T).5*dx,ymin-(T).5*dy,ymax+(T).5*dy,true);}

    GRID_1D<T> Remove_Dimension(int dimension) const
    {return GRID_1D<T>(Counts().Remove_Index(dimension),Domain().Remove_Dimension(dimension),Is_MAC_Grid());}

    GRID_1D<T> Get_Horizontal_Grid() const
    {return GRID_1D<T>(m,xmin,xmax,Is_MAC_Grid());}

    GRID_1D<T> Get_1D_Grid(const int axis) const
    {return GRID_1D<T>(Counts()[axis],Domain().Minimum_Corner()[axis],Domain().Maximum_Corner()[axis],Is_MAC_Grid());}

    static TV_INT Node_Neighbor(const TV_INT& index,const int i)
    {static const TV_INT neighbor_offset[4]={TV_INT(-1,0),TV_INT(1,0),TV_INT(0,-1),TV_INT(0,1)};
    assert(1<=i&&i<=4);return index+neighbor_offset[i-1];}

    void Cells_Neighboring_Node(const TV_INT& node,TV_INT cells[4]) const
    {cells[0]=TV_INT(node.x-1,node.y-1);cells[1]=TV_INT(node.x,node.y-1);cells[2]=TV_INT(node.x-1,node.y);cells[3]=TV_INT(node.x,node.y);}

    void Nodes_In_Cell_From_Minimum_Corner_Node(const TV_INT& minimum_corner_node,TV_INT nodes[4]) const
    {nodes[0]=minimum_corner_node;nodes[1]=TV_INT(minimum_corner_node.x+1,minimum_corner_node.y);nodes[2]=TV_INT(minimum_corner_node.x,minimum_corner_node.y+1);
    nodes[3]=TV_INT(minimum_corner_node.x+1,minimum_corner_node.y+1);}

    void Node_Locations_In_Cell_From_Minimum_Corner_Node(const TV_INT& minimum_corner_node,ARRAY<TV>& node_locations) const
    {assert(node_locations.m==number_of_nodes_per_cell);
    TV_INT nodes[number_of_nodes_per_cell];Nodes_In_Cell_From_Minimum_Corner_Node(minimum_corner_node,nodes);
    for(int i=1;i<=number_of_nodes_per_cell;i++) node_locations(i)=Node(nodes[i-1]);}

    TV_INT First_Face_Index_In_Cell(const int axis,const TV_INT& cell_index) const
    {return cell_index;}

    TV_INT Second_Face_Index_In_Cell(const int axis,const TV_INT& cell_index) const
    {TV_INT face_index=cell_index;face_index[axis]++;return face_index;}

    static TV_INT One_Ring_Neighbor(const TV_INT& index,const int i)
    {static const TV_INT neighbor_offset[8]={TV_INT(-1,-1),TV_INT(0,-1),TV_INT(1,-1),TV_INT(-1,0),TV_INT(1,0),TV_INT(-1,1),
        TV_INT(0,1),TV_INT(1,1)};
    assert(1<=i&&i<=8);return index+neighbor_offset[i-1];}

//#####################################################################
    void Initialize(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid=false);
    static GRID_2D<T> Create_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes=0);
    static GRID_2D<T> Create_Even_Sized_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes=0);
    template<class RW> void Read(std::istream& input);
    template<class RW> void Write(std::ostream& output) const;
//#####################################################################
};
// global functions
template<class T>
inline std::ostream& operator<<(std::ostream& output,const GRID_2D<T>& grid)
{output << grid.Domain() << " divided by ("<<grid.m<<","<<grid.n<<")";return output;}
}
#include <Arrays/ARRAYS_2D.h>
#endif
