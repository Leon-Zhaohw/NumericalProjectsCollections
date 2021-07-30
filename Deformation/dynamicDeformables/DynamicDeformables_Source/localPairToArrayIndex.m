uchar FixedElement::localPairToArrayIndex(uchar i, uchar j)
{
    // assert is here for the cautious. 
    // If you pass negative indices that's 
    // another thing you should never do! 
    assert(i <= j);
    uchar n = (_numParticlesInElement << 1) - 1;
    return j + ((n - i ) * i) >> 1;
}