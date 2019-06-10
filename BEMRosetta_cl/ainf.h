//// Infinite added-mass (Ainf) calculator
//
//   This function calculates the added-mass radiation coefficient at the
//   infinite-frequency (Ainf) using the Ogilvie's formula [A], which is 
//   necessary for time-domain simulations of semi-submerged or totally submerged bodies.
//   
//   	    [A] Ogilvie, F. T. Recent Progress Towards the Understanding and Prediction of Ship
//   		5th Symposium on Naval Hydrodynamics, vol. 112, Washington DC, USA, 1964.
//
//   Inputs: w --> vector of frequencies [rad/s]     --> n_wx1
//           A --> Radiation added-mass vector [kg]  --> n_wx1
//           B --> Radiation damping vector [Ns/m]   --> n_wx1
//           T --> Time-window for the calculation of the radiation impulse response function --> n_tx1
//		where n_w and n_t are the number of frequencies and length of the time vector, respectively.
//
//   Outputs: Ainf --> radiation added-mass at infinite frequency
//
// Markel Penalba
// Fluid-Mechanics research group, Mondragon University

template <class Range, class V>
V Ainf(const Range &w, const Range &A, const Range &B, const Range &T) {
    int n = T.GetCount();   	// Length of the time-vector
    int nw = w.GetCount();   	// Length of the vector of frequencies
    
    Buffer<V> K(n, 0); 			// Allocation of the impulse response vector
    for (int ki = 0; ki < n; ++ki) {
        for (int kj = 1; kj < nw; ++kj)
            K(ki) += B(kj)*cos(w[kj]*T[ki])*(w[kj] - w[kj-1])*2/M_PI ;
    }
    
    Buffer<V> Kint(nw, 0);		// Allocation of the vector with the impulse response integral
    V Ainf = 0;
    for (int kj = 0; kj < nw; ++kj) {
        for (int ki = 1; ki < n; ++ki) 
            Kint[kj] += K[ki]*sin(w[kj]*T[ki])*(T[ki] - T[ki-1]);
        // Ogilvie's formula
        Kint[kj] = A[kj] + Kint[kj]/w[kj];
        Ainf += Kint[kj];
	}
    return Ainf/nw;
}