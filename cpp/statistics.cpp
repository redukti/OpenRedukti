/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#include <statistics.h>

#include <algorithm>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace redukti
{

double mean(const double *start, const double *end)
{
	assert(start < end);

	double m = 0.0;
	int k = 0;
	for (const double *p = start; p < end; p++, k++) {
		m += *p;
	}
	return m / k;
}

double variance(const double *start, const double *end, bool unbiased)
{
	assert(start < end);

	double s = 0.0;
	double m_kminus1 = *start;
	int k = 1;
	for (const double *p = start + 1; p < end; p++) {
		k += 1;
		double m_k = m_kminus1 + (*p - m_kminus1) / k;
		s = s + (*p - m_kminus1) * (*p - m_k);
		m_kminus1 = m_k;
	}
	return s / (k - (unbiased ? 1 : 0));
}

double standard_deviation(const double *start, const double *end, bool unbiased)
{
	return sqrt(variance(start, end, unbiased));
}

struct abs_diff {
	double operator()(double a, double b) const { return b - a; }
};
struct rel_diff {
	double operator()(double a, double b) const { return (b - a) / a; }
};

template <typename diff_type>
void calculate_difference(const double *inp_start, const double *inp_end, double *outp_iter, int lag,
			  const diff_type &diff_op)
{
	for (const double *p = inp_start + lag; p < inp_end; p++, inp_start++, outp_iter++)
		*outp_iter = diff_op(*inp_start, *p);
}

void lagged_difference(double *dest, const double *src, size_t src_size, int lag, DiffType difference_type)
{
	assert(lag > 0 && lag < src_size);
	if (difference_type == DIFFERENCE_ABSOLUTE)
		return calculate_difference(src, src + src_size, dest, lag, abs_diff());
	else
		return calculate_difference(src, src + src_size, dest, lag, rel_diff());
}

struct sd_func {
	double operator()(double prev, double v, double lambda) const
	{
		return sqrt(lambda * prev * prev + v * v * (1.0 - lambda));
	}
};

struct var_func {
	double operator()(double prev, double v, double lambda) const { return lambda * prev + v * v * (1.0 - lambda); }
};

struct abs_func {
	double operator()(double prev, double v, double lambda) const
	{
		return lambda * prev + fabs(v) * (1.0 - lambda);
	}
};

template <typename stat_func>
static void ewma(const double *start, const double *end, double *dst, double lambda, double seed, const stat_func &f)
{
	for (const double *p = start; p != end; p++, dst++) {
		seed = f(seed, *p, lambda);
		*dst = seed;
	}
}

void ewma_volatility(double *dst, const double *start, size_t src_size, double lambda, double seed,
		     VolatilityEstimator method)
{
	if (method == STANDARD_DEVIATION)
		ewma(start, start + src_size, dst, lambda, seed, sd_func());
	else if (method == VARIANCE)
		ewma(start, start + src_size, dst, lambda, seed, var_func());
	else
		ewma(start, start + src_size, dst, lambda, seed, abs_func());
}

int equal(double d1, double d2) { return fabs(d1 - d2) < 1e-6; }
int equal3(double d1, double d2) { return fabs(d1 - d2) < 1e-3; }

static int test_basics()
{
	double data[] = {3.6, 1.8, 3.333, 2.283, 4.533, 2.883};
	double v = variance(std::begin(data), std::end(data), true);
	if (!equal(v, 0.951530))
		return 1;
	double v1[5], v2[4];
	v1[0] = 3;
	v1[1] = 5;
	v1[2] = 8;
	v1[3] = 12;
	v1[4] = 17;
	lagged_difference(v2, v1, 5, 1);
	if (v2[0] != v1[1] - v1[0] || v2[1] != v1[2] - v1[1] || v2[2] != v1[3] - v1[2] || v2[3] != v1[4] - v1[3])
		return 1;
	return 0;
}

static int test_ewma()
{
	double data[9]{22.0, 22.5, 21.0, 23.0, 25, 23.5, 23, 22, 23.6};
	double diffdata[8];
	double expecteddif[8] = {0.5, -1.5, 2.0, 2.0, -1.5, -0.5, -1.0, 1.6};
	lagged_difference(diffdata, data, 9, 1, DIFFERENCE_ABSOLUTE);
	for (int i = 0; i < 8; i++) {
		if (!equal(diffdata[i], expecteddif[i]))
			return 1;
	}

	double m = mean(std::begin(diffdata), std::end(diffdata));
	if (!equal(m, 0.2))
		return 1;

	double seed = standard_deviation(std::begin(diffdata), std::end(diffdata), true);
	if (!equal(seed, 1.523154621))
		return 1;

	double ewma[8] = {0};
	double ewma_expect[8] = {1.502631, 1.502552, 1.519846, 1.536436, 1.535355, 1.514628, 1.501757, 1.504797};
	ewma_volatility(ewma, std::begin(diffdata), 8, 0.97, seed, STANDARD_DEVIATION);
	for (int i = 0; i < 8; i++) {
		if (!equal(ewma_expect[i], ewma[i]))
			return 1;
	}

	// This test is taken from example in Risk Metrics Tech document
	double rmeg[20] = {0.633,  0.115,  -0.459, 0.093, 0.176, -0.087, -0.142, 0.324,  -0.943, -0.528,
			   -0.107, -0.159, -0.445, 0.053, 0.152, -0.318, 0.424,  -0.708, -0.105, -0.257};
	double expected[20] = {0.401, 0.378, 0.368, 0.346, 0.327, 0.308, 0.291, 0.280, 0.316, 0.314,
			       0.296, 0.280, 0.275, 0.258, 0.244, 0.236, 0.232, 0.248, 0.234, 0.224};
	double got[20] = {0};
	seed = rmeg[0] * rmeg[0];
	ewma_volatility(got, std::begin(rmeg), 20, 0.94, seed, VARIANCE);
	for (int i = 0; i < 20; i++) {
		if (!equal3(expected[i], got[i]))
			return 1;
	}
	return 0;
}

int test_statistics()
{
	int failure_count = 0;
	failure_count += test_basics();
	failure_count += test_ewma();
	if (failure_count != 0)
		printf("Statistics Tests FAILED\n");
	else
		printf("Statistics Tests OK\n");
	return failure_count;
}

} // namespace redukti
