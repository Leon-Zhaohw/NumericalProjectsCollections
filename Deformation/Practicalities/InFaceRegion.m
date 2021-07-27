bool InFaceRegion(p, triFace)
   
    // If the projection of the point into the face plane 
    // resides inside the face boundary 
    if pproj in triFace:
        return true
    
    distToPlane[3] = {triFace.DistToBisector(0,p),
                      triFace.DistToBisector(1,p),
                      triFace.DistToBisector(2,p)}
    
    // Test if the point is in the cell of the face
    if all distToPlane >= 0:
        return true
    
    // Is the point outside  the plane of a boundary edge, 
    // and inside the planes of all non-boundary edges?   
    if (triFace.IsBoundaryEdge(0) and distToPlane[0] < 0) or 
       (triFace.IsBoundaryEdge(1) and distToPlane[1] < 0) or  
       (triFace.IsBoundaryEdge(2) and distToPlane[2] < 0):
       
        for e in [0,2]:
            if not triFace.IsBoudaryEdge(e) and distToPlane[e] < 0:
                return false
	
        return true
    
  // The point is not in the region of the face
  return false