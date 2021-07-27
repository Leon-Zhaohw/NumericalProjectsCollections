bool HandleContinuousCollisions(CollisionDetector, maxIters)
    
    CollisionDetector->DetectInitialCollisions()
    collisions = CollisionDetector->GetResults()    
    iter = 0
    bool collisionFree = true
    bool allRegionsResolved = false
    if not collisions.empty():
        collisionFree = false
        GatherCollisionsPerDynamicVertex(collisions)
        CreateContactRegions(collisions)
        
        while !collisionFree and iter < maxIters: 
            // Solve first contacts in each disjoint regions
            // and update the end position state marking 
            // incident mesh features for the BVH update.
            SolveContactRegions()            
            iters++
            allRegionsResolved = AreAllContactRegionsCollisionFree()
            // Steps for incremental detection 
            CollisionDetector->UpdateBVHsWithMarkedFeatures()
            CollisionDetector->PerformIncrementalDetection()            
            CollisionDetector->GetResults(newCollisions)

            // Determine whether collisions in newCollisions
            // are genuinely new or are already represented 
            // in the contact regions. Remove already discovered
            // collisions from newCollisions.
            ProcessContacts(newCollisions)
            if not newCollisions.empty():         
                // Contact regions may just expand
                // or expand and merge with others   
                ExpandAndMergeContactRegions(newCollisions)             
            else if allRegionsResolved:              
                collisionFree = true           

    return collisionFree
