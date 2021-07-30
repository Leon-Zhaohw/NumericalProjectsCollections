```
void FLUID_3D_DCT::projectDCT() {
  Timer timer;
  timer.Reset();
  
  ForwardTransformToFrequency(_velocity, vxFreq_, vyFreq_, vzFreq_);
  totalDCTtime_ += timer.ElapsedTimeInSeconds();
  
  // Project the component that along the wavenumber away.
  for (int kz = 0; kz < _zRes; kz++) {
    for (int ky = 0; ky < _yRes; ky++) {
      for (int kx = 0; kx < _xRes; kx++) {
        const int index = kx + ky*_xRes + kz*_slabSize;
        
        // Wavenumber.
        Eigen::Vector3d K(static_cast<double>(kx)*invxRes_,
                          static_cast<double>(ky)*invyRes_,
                          static_cast<double>(kz)*invzRes_);
        Eigen::Vector3d W; // (vxFreq_[index], vyFreq_[index], vzFreq_[index]);
        if (kx == 0) W[0] = 0; else W[0] = vxFreq_[index - 1];
        if (ky == 0) W[1] = 0; else W[1] = vyFreq_[index - _xRes];
        if (kz == 0) W[2] = 0; else W[2] = vzFreq_[index - _slabSize];
        
        if (kz != 0 || ky!= 0 || kx != 0) {
          W = W - 1.0 / K.squaredNorm() * (K.dot(W))*K;
        }
        
        if (kx != 0) vxFreq_[index - 1] = W[0];
        if (ky != 0) vyFreq_[index - _xRes] = W[1];
        if (kz != 0) vzFreq_[index - _slabSize] = W[2];
      }
    }
  }
  // Inverse transform to velocity.
  timer.Reset();
  InverseTransformToVelocity(vxFreq_, vyFreq_, vzFreq_, &_velocity);
  totalDCTtime_ += timer.ElapsedTimeInSeconds();
}
```

