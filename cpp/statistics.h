/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017-2019 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#ifndef _REDUKTI_STATISTICS_H
#define _REDUKTI_STATISTICS_H

#include <stdio.h>

namespace redukti
{

extern double mean(const double *start, const double *end);
extern double variance(const double *start, const double *end, bool unbiased = true);
extern double standard_deviation(const double *start, const double *end, bool unbiased = true);

// For every column in m, perform a diff between element at i+lag and i,
// an return a new matrix.

enum DiffType {
	DIFFERENCE_ABSOLUTE = 1, // Absolute Differences
	DIFFERENCE_RELATIVE = 2  // Relative Differences
};

// The src_len parameter should say how big src is
// The dest is assumed to be src_len-lag in size
extern void lagged_difference(double *dest, const double *src, size_t src_len, int lag,
			      DiffType difference_type = DIFFERENCE_ABSOLUTE);

enum VolatilityEstimator {
	MEAN_ABSOLUTE_DEVIATION = 1, // Mean Absolute Deviation
	STANDARD_DEVIATION = 2,      // Standard Deviation
	VARIANCE = 3		     // Variance
};

// Both dst and src must be of len n
// The seed must be compatible with the type of dispersion being requested
extern void ewma_volatility(double *dst, const double *src, size_t n, double lambda, double seed,
			    VolatilityEstimator method = STANDARD_DEVIATION);

extern int test_statistics();

} // namespace redukti

#endif
