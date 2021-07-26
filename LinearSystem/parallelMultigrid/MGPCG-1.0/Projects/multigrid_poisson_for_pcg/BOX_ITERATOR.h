//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Yongning Zhu, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
// Class BOX_ITERATOR
//#####################################################################
#ifndef __BOX_ITERATOR__
#define __BOX_ITERATOR__

namespace PhysBAM{

template<int d,int stride=1> class BOX_ITERATOR;

template<int d> 
class BOX_ITERATOR<d,1>
{
    typedef VECTOR<int,d> TV_INT;
    typedef BOX<TV_INT> T_BOX;

    const T_BOX& box;
    TV_INT index;

public:
    BOX_ITERATOR(const T_BOX& box_input)
        :box(box_input)
    {
        Reset();
    }

    void Reset()
    {index=box.min_corner;}

    bool Valid() const
    {return index.x<=box.max_corner.x;}

    void Next()
    {for(int i=d;i>=1;i--) if(index(i)<box.max_corner(i) || i==1){index(i)++;return;} else index(i)=box.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

template<int d,int stride>
class BOX_ITERATOR
{
    STATIC_ASSERT((stride!=1));
    typedef VECTOR<int,d> TV_INT;
    typedef BOX<TV_INT> T_BOX;

    const T_BOX& box;
    TV_INT index;

public:
    BOX_ITERATOR(const T_BOX& box_input)
        :box(box_input)
    {
        Reset();
    }

    void Reset()
    {index=box.min_corner;}

    bool Valid() const
    {return index.x<=box.max_corner.x;}

    void Next()
    {for(int i=d;i>=1;i--) if(index(i)+stride<=box.max_corner(i) || i==1){index(i)+=stride;return;} else index(i)=box.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

template<int d>
class BOUNDARY_ITERATOR
{
    typedef VECTOR<int,d> TV_INT;
    typedef BOX<TV_INT> T_BOX;
    
    const T_BOX& box;
    TV_INT index;
    T_BOX regions[2*d];

    int current_region;
    int side; // side=1: xmin, side=2 xmax, side=3:ymin, side=4:ymax
    bool valid;
public:
    BOUNDARY_ITERATOR(const T_BOX& box_input,int side=0)
	:box(box_input),side(side)
    {
	assert(side<=2*d);
	TV_INT min_corner=box.min_corner+1;
	TV_INT max_corner=box.max_corner-1;
	for(int v=1;v<=d;v++){
	    min_corner(v)=box.min_corner(v);
	    max_corner(v)=box.min_corner(v);
	    regions[(v-1)*2]=T_BOX(min_corner,max_corner);
	    min_corner(v)=box.max_corner(v);
	    max_corner(v)=box.max_corner(v);
	    regions[(v-1)*2+1]=T_BOX(min_corner,max_corner);
	    min_corner(v)=box.min_corner(v);
	}
	if(side)
	    Reset(side-1);
	else
	    Reset();
    }	
    
    void Reset(const int region_index=0)
    {
	valid=side?region_index<side:(region_index<2*d);
	if(valid){
	    current_region=region_index;
	    index=regions[current_region].min_corner;
	}
    }

    bool Valid() const
    {
	return valid;
    }

    void Next()
    {
	if(!valid) return;
	for(int i=d;i>=1;i--){
	    if(index(i)<regions[current_region].max_corner(i)){
		index(i)++; return;
	    }
	    else
		index(i)=regions[current_region].min_corner(i);
	}
	Reset(current_region+1);
    }
	
    const TV_INT&  Index()
    {return index;}
	

};

//#####################################################################
}
#endif
