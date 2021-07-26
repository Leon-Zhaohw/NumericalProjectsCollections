//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Duc Nguyen, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_1D
//#####################################################################
// A number of functions (e.g. the Clamp functions) assume the grid indexing starts at 1.  This way we can use truncation rather than floor because floor is really slow.
//#####################################################################
#ifndef __GRID_1D__
#define __GRID_1D__

#include <Grids/POLICY_UNIFORM.h>
#include <Arrays/ARRAYS_FORWARD.h>
#include <Geometry/BOX.h>
#include <Grids/GRID_0D.h>
namespace PhysBAM{

template<class T> class UNIFORM_GRID_ITERATOR_NODE_1D;
template<class T> class UNIFORM_GRID_ITERATOR_CELL_1D;
template<class T> class UNIFORM_GRID_ITERATOR_FACE_1D;
template<class T> class LEVELSET_1D;
template<class T> class INCOMPRESSIBLE_1D;
template<class T_GRID> class LINEAR_INTERPOLATION_MAC_1D_HELPER;
class FLOOD_FILL_1D;

template<class T>
class GRID_1D:public POLICY_UNIFORM<VECTOR<T,1> >
{
    STATIC_ASSERT((IS_SCALAR<T>::value));
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
public:
    enum REGION {WHOLE_REGION,GHOST_REGION,BOUNDARY_REGION,INTERIOR_REGION,BOUNDARY_INTERIOR_REGION}; // for iterators
    static const int dimension=1;
    static const int number_of_cells_per_block=2;
    static const int number_of_faces_per_block=3;
    static const int number_of_incident_faces_per_block=1;
    static const int number_of_nodes_per_face=1;
    static const int number_of_nodes_per_cell=2;
    static const int number_of_cells_per_node=2;
    static const int number_of_neighbors_per_node=2;
    static const int number_of_neighbors_per_cell=2;
    static const int number_of_faces_per_cell=2;
    static const int number_of_one_ring_neighbors_per_cell=2;
    typedef UNIFORM_GRID_ITERATOR_NODE_1D<T> NODE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_CELL_1D<T> CELL_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_FACE_1D<T> FACE_ITERATOR;
    typedef LEVELSET_1D<T> LEVELSET;
    typedef INCOMPRESSIBLE_1D<T> INCOMPRESSIBLE;
    typedef FLOOD_FILL_1D FLOOD_FILL;
    typedef LINEAR_INTERPOLATION_MAC_1D_HELPER<GRID_1D> LINEAR_INTERPOLATION_MAC_HELPER;

    int m; // # of points: x direction
    T xmin,xmax; // left and right wall
    T dx; // cell size
    T one_over_dx;
    T MAC_offset; // 0 for a regular grid and .5 for a MAC grid
    int number_of_cells_x; // not saved to file

    GRID_1D()
    {
        Initialize(TV_INT(),BOX<TV>::Unit_Box());
    }

    GRID_1D(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid=false)
    {
        Initialize(m_input,xmin_input,xmax_input,MAC_grid);
    }

    GRID_1D(const int m_input,const BOX<TV>& box,const bool MAC_grid=false)
    {
        Initialize(m_input,box.min_corner.x,box.max_corner.x,MAC_grid);
    }

    GRID_1D(const TV_INT& counts,const BOX<TV>& box,const bool MAC_grid=false)
    {
        Initialize(counts,box,MAC_grid);
    }

    GRID_1D(const GRID_1D<T>& grid_input)
    {
        Initialize(grid_input.Counts(),grid_input.Domain(),grid_input.Is_MAC_Grid());
    }

    template<class T2>
    GRID_1D(const GRID_1D<T2>& grid_input)
    {
        Initialize(grid_input.Counts(),BOX<TV>(grid_input.Domain()),grid_input.Is_MAC_Grid());
    }

    GRID_1D<T>& operator=(const GRID_1D<T>& grid_input)
    {if(this!=&grid_input) Initialize(grid_input.Counts(),grid_input.Domain(),grid_input.Is_MAC_Grid());return *this;}

    void Set_To_Double_Resolution_Grid(const GRID_1D<T>& grid)
    {Initialize(2*grid.Numbers_Of_Cells()+!grid.Is_MAC_Grid(),grid.Domain(),grid.Is_MAC_Grid());}

    void Initialize(const TV_INT& counts,const BOX<TV>& box,const bool MAC_grid=false)
    {Initialize(counts.x,box,MAC_grid);}

    void Initialize(const int m_input,const BOX<TV>& box,const bool MAC_grid=false)
    {Initialize(m_input,box.min_corner.x,box.max_corner.x,MAC_grid);}

    bool Is_MAC_Grid() const
    {return MAC_offset==.5;}

    TV DX() const
    {return TV(dx);}

    TV One_Over_DX() const
    {return TV(one_over_dx);}

    TV_INT Numbers_Of_Cells() const
    {return TV_INT(number_of_cells_x);}

    TV_INT Numbers_Of_Nodes() const
    {return Numbers_Of_Cells()+1;}

    T Minimum_Edge_Length() const
    {return dx;}

    T Cell_Diagonal() const
    {return dx;}

    T Maximum_Edge_Length() const
    {return dx;}

    T Cell_Size() const
    {return dx;}

    T Face_Size(const int axis) const
    {assert(axis==1);return 1;}

    TV Xmin() const
    {return TV(xmin);}

    TV Xmax() const
    {return TV(xmax);}

    T x(const int i) const // x at the i=1 to i=m grid points
    {return xmin+((T)i-(T)1+MAC_offset)*dx;}

    T x_plus_half(const int i) const
    {return xmin+((T)i-(T).5+MAC_offset)*dx;}

    T x_minus_half(const int i) const
    {return xmin+((T)i-(T)1.5+MAC_offset)*dx;}

    TV X(const int i) const
    {return TV(x(i));}

    TV X(const TV_INT& index) const
    {return TV(x(index.x));}

    TV Node(const int i) const
    {return TV(xmin+(i-(T)1)*dx);}

    TV Node(const TV_INT& index) const
    {return TV(xmin+(index.x-(T)1)*dx);}

    TV Center(const int i) const
    {return TV(xmin+(i-(T).5)*dx);}

    TV Center(const TV_INT& index) const
    {return TV(xmin+(index.x-(T).5)*dx);}

    TV X_Face(const int i) const
    {return TV(xmin+(i-(T)1)*dx);}
    
    TV X_Face(const TV_INT& index) const
    {return TV(xmin+(index.x-(T)1)*dx);}

    TV Face(const int axis,const TV_INT& index) const
    {assert(axis==1);return TV(xmin+(index.x-(T)1)*dx);}

    BOX<TV> Face_Domain(const int axis,const TV_INT& index,const T thickness_over_two=0) const
    {TV dimensions=(T).5*DX();dimensions[axis]=thickness_over_two;BOX<TV> domain(Face(axis,index));domain.Change_Size(dimensions);return domain;}

    static TV_INT Node_Cell_Index(const TV_INT& node,const int cell)
    {static const TV_INT cell_from_node_offset[2]={TV_INT(-1),TV_INT(0)};
    assert(1<=cell&&cell<=2);return node+cell_from_node_offset[cell-1];}

    static TV_INT Face_Node_Index(const int axis,const TV_INT& face_index,const int node)
    {assert(axis==1&&node==1);return face_index;}

    static TV_INT Node_Face_Index(const int axis,const TV_INT& node_index,const int face)
    {assert(axis==1&&face==1);return node_index;}

    void Index(const TV& location,int& i) const // returns the left index
    {i=(int)floor((location.x-xmin)*one_over_dx+1-MAC_offset);} // note that "floor" is expensive

    void Cell(const TV& location,int& i,const int number_of_ghost_cells) const // returns the left
    {i=(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells+1)-number_of_ghost_cells;}

    TV_INT Cell(const TV& location,const int number_of_ghost_cells) const // returns the left
    {TV_INT index;Cell(location,index.x,number_of_ghost_cells);return index;}

    BOX<TV> Cell_Domain(const int i) const
    {return BOX<TV>(xmin+(i-(T)1)*dx,xmin+i*dx);}

    BOX<TV> Cell_Domain(const TV_INT& index) const
    {return BOX<TV>(xmin+(index.x-(T)1)*dx,xmin+index.x*dx);}

    static void Cells_Touching_Face(const int axis,const TV_INT& face_index,TV_INT& cell1,TV_INT& cell2)
    {cell2=face_index;cell1.x=face_index.x-1;}

    TV Clamp(const TV& location) const
    {return TV(clamp(location.x,xmin,xmax));}

    void Clamp(int& i) const
    {i=clamp(i,1,m);}

    void Clamp(TV_INT& index) const
    {index.x=clamp(index.x,1,m);}

    void Clamp_End_Minus_One(int& i) const
    {i=clamp(i,1,m-1);}

    void Clamped_Index(const TV& location,int& i) const
    {i=min(m,1+max(0,(int)((location.x-xmin)*one_over_dx-MAC_offset)));}

    TV_INT Clamped_Index(const TV& location) const
    {return TV_INT(min(m,1+max(0,(int)((location.x-xmin)*one_over_dx-MAC_offset))));}

    TV_INT Clamped_Index_End_Minus_One(const TV& location) const
    {return TV_INT(min(m-1,1+max(0,(int)((location.x-xmin)*one_over_dx-MAC_offset))));}

    void Clamped_Index_End_Minus_One(const TV& location,int& i,const int number_of_ghost_cells) const
    {i=min(m+number_of_ghost_cells-1,1-number_of_ghost_cells+max(0,(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells-MAC_offset)));}

    TV_INT Clamped_Index_End_Minus_One(const TV& location,const int number_of_ghost_cells) const
    {TV_INT index;Clamped_Index_End_Minus_One(location,index.x,number_of_ghost_cells);return index;}

    void Clamp_To_Cell(const TV& location,int& i) const
    {i=min(number_of_cells_x,1+max(0,(int)((location.x-xmin)*one_over_dx)));}

    TV_INT Clamp_To_Cell(const TV& location) const
    {return TV_INT(min(number_of_cells_x,1+max(0,(int)((location.x-xmin)*one_over_dx))));}

    void Clamp_To_Cell(const TV& location,int& i,const int number_of_ghost_cells) const
    {i=min(number_of_cells_x+number_of_ghost_cells,1-number_of_ghost_cells+max(0,(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells)));}

    TV_INT Clamp_To_Cell(const TV& location,const int number_of_ghost_cells) const
    {return TV_INT(min(number_of_cells_x+number_of_ghost_cells,1-number_of_ghost_cells+max(0,(int)((location.x-xmin)*one_over_dx+number_of_ghost_cells))));}

    BOX<TV_INT> Clamp_To_Cell(const BOX<TV>& box,const int number_of_ghost_cells) const
    {return BOX<TV_INT>(Clamp_To_Cell(box.Minimum_Corner(),number_of_ghost_cells),Clamp_To_Cell(box.Maximum_Corner(),number_of_ghost_cells));}

    TV_INT Block_Index(const TV& X,const int number_of_ghost_cells) const // index of node at center of block
    {assert(Is_MAC_Grid());return Clamped_Index_End_Minus_One(X,number_of_ghost_cells)+TV_INT::All_Ones_Vector();}

    void Closest_Node(const TV& location,int& i) const
    {i=min(m,1+max(0,(int)((location.x-xmin)*one_over_dx+(T).5)));}

    TV_INT Closest_Node(const TV& location) const
    {TV_INT index;Closest_Node(location,index.x);return index;}

    BOX<TV> Domain() const
    {return BOX<TV>(xmin,xmax);}

    BOX<TV> Ghost_Domain(const int number_of_ghost_cells) const
    {return BOX<TV>(xmin-dx*number_of_ghost_cells,xmax+dx*number_of_ghost_cells);}

    TV_INT Counts() const
    {return TV_INT(m);}

    BOX<TV_INT> Domain_Indices() const
    {return BOX<TV_INT>(1,m);}

    BOX<TV_INT> Node_Indices(const int ghost_cells=0) const
    {return BOX<TV_INT>(TV_INT()+1,Numbers_Of_Nodes()).Thickened(ghost_cells);}

    BOX<TV_INT> Cell_Indices(const int ghost_cells=0) const
    {return BOX<TV_INT>(TV_INT()+1,Numbers_Of_Cells()).Thickened(ghost_cells);}

    BOX<TV_INT> Block_Indices(const int ghost_cells=0) const
    {return Node_Indices(ghost_cells);}

    VECTOR<BOX<TV_INT>,1> Face_Indices(const int ghost_cells=0) const
    {return VECTOR<BOX<TV_INT>,1>(Get_X_Face_Grid().Node_Indices(ghost_cells));}

    bool Outside(const T location) const
    {return location < xmin || location > xmax;}

    bool Outside(const TV& location) const
    {return location.x < xmin || location.x > xmax;}

    GRID_1D<T> Get_MAC_Grid() const
    {return GRID_1D<T>(Numbers_Of_Cells(),Domain(),true);}

    GRID_1D<T> Get_Regular_Grid() const
    {return GRID_1D<T>(Numbers_Of_Nodes(),Domain(),false);}

    GRID_1D<T> Get_X_Face_Grid() const
    {return GRID_1D<T>(number_of_cells_x+1,xmin,xmax);}

    GRID_1D<T> Get_Face_Grid(const int axis) const
    {assert(axis==1);return Get_X_Face_Grid();}

    GRID_1D<T> Get_Regular_Grid_At_MAC_Positions() const
    {assert(Is_MAC_Grid());
    return GRID_1D<T>(m,xmin+(T).5*dx,xmax-(T).5*dx);}

    GRID_1D<T> Get_MAC_Grid_At_Regular_Positions() const
    {assert(!Is_MAC_Grid());return GRID_1D<T>(m,xmin-(T).5*dx,xmax+(T).5*dx,true);}

    GRID_0D<T> Remove_Dimension(int dimension) const
    {return GRID_0D<T>(Counts().Remove_Index(dimension),Domain().Remove_Dimension(dimension),Is_MAC_Grid());}

    GRID_1D<T> Get_1D_Grid(int axis) const
    {return *this;}

    static TV_INT Node_Neighbor(const TV_INT& index,const int i) // i=1 to 2
    {static const TV_INT neighbor_offset[2]={TV_INT(-1),TV_INT(1)};
    assert(1<=i&&i<=2);return index+neighbor_offset[i-1];}

    void Cells_Neighboring_Node(const TV_INT& node,TV_INT cells[2]) const
    {cells[0]=TV_INT(node.x-1);cells[1]=TV_INT(node.x);}

    void Nodes_In_Cell_From_Minimum_Corner_Node(const TV_INT& minimum_corner_node,TV_INT nodes[2]) const
    {nodes[0]=minimum_corner_node;nodes[1]=TV_INT(minimum_corner_node.x+1);}

    void Node_Locations_In_Cell_From_Minimum_Corner_Node(const TV_INT& minimum_corner_node,ARRAY<TV>& node_locations) const
    {assert(node_locations.m==number_of_nodes_per_cell);
    TV_INT nodes[number_of_nodes_per_cell];Nodes_In_Cell_From_Minimum_Corner_Node(minimum_corner_node,nodes);
    for(int i=1;i<=number_of_nodes_per_cell;i++) node_locations(i)=Node(nodes[i-1]);}

    TV_INT First_Face_Index_In_Cell(const int axis,const TV_INT& cell_index) const
    {return cell_index;}

    TV_INT Second_Face_Index_In_Cell(const int axis,const TV_INT& cell_index) const
    {TV_INT face_index=cell_index;face_index[axis]++;return face_index;}

    static TV_INT One_Ring_Neighbor(const TV_INT& index,const int i)
    {return Node_Neighbor(index,i);}

//#####################################################################
    void Initialize(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid=false);
    static GRID_1D<T> Create_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes=0);
    static GRID_1D<T> Create_Even_Sized_Grid_Given_Cell_Size(const BOX<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes=0);
    template<class RW> void Read(std::istream& input);
    template<class RW> void Write(std::ostream& output) const;
//#####################################################################
};
// global functions
template<class T>
inline std::ostream& operator<<(std::ostream& output,const GRID_1D<T>& grid)
{output << grid.Domain() << " divided by ("<<grid.m<<")";return output;}
}
#include <Arrays/ARRAY.h>
#include <Arrays/ARRAYS_1D.h>
#endif
