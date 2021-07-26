/*
 *  sampleseeder.h
 *  eltopo2d_project
 *
 *  Created by tyson on 17/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "vec.h"

namespace eltopo2d
{
   class DynamicSurface;
}


class SampleSeeder
{
   
public:
   
   static void generate_samples( const eltopo2d::DynamicSurface& surface, 
                                 double desired_dx, 
                                 std::vector<Vec2f>& samples,
                                 std::vector<int>& expected_signs);
   
   static void generate_regular_samples( const Vec2d& domain_min,
                                         const Vec2d& domain_max,
                                         double dx,
                                         std::vector<Vec2f>& samples,
                                         std::vector<int>& expected_signs );
   
};