rm -f meshes.h
../../skinning/fixObj --dump-arrays univJointLeft univJointLeft.obj >> meshes.h
../../skinning/fixObj --dump-arrays univJointRight univJointRight.obj >> meshes.h
../../skinning/fixObj --dump-arrays hingeJointLeft hingeJointLeft.obj >> meshes.h
../../skinning/fixObj --dump-arrays hingeJointRight hingeJointRight.obj >> meshes.h
../../skinning/fixObj --dump-arrays ballJointLeft ballJointLeft.obj >> meshes.h
../../skinning/fixObj --dump-arrays ballJointRight ballJointRight.obj >> meshes.h

