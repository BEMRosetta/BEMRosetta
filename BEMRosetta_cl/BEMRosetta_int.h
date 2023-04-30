// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 - 2022, the BEMRosetta author and contributors
#ifndef _BEMRosetta_cl_BEMRosetta_int_h_
#define _BEMRosetta_cl_BEMRosetta_int_h_

#include "functions.h"

int IsTabSpace(int c);

template <class Range>
void AvgB(Range &ret, const UArray<const Range*> &d) {
	int numT = d.size();
	if (numT == 0) 
		return;
	
	int num = int(d[0]->size());
	for (int it = 0; it < numT; ++it) 
		if (d[it]->size() != num)
			throw Exc(t_("Avg() has to have same number of values"));
	
	for (int i = 0; i < num; ++i) {
		Eigen::VectorXd r(numT);
		for (int it = 0; it < numT; ++it) 
			r[it] = (*d[it])[i];
		ret[i] = r.mean();
	}
}


#endif
