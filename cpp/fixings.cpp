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

#include <dayfractions.h>
#include <fixings.h>
#include <index.h>
#include <timeseries.h>

namespace redukti
{

static int test_timeseries()
{
	Value values[] = {
	    Value(make_date(1, January, 2014), 5.0),
	    Value(make_date(30, December, 2013), 4.9),
	    Value(make_date(2, January, 2014), 5.1),
	    Value(make_date(3, January, 2014), 5.2),
	};
	Value additional_values[] = {
	    Value(make_date(29, December, 2013), 4.81),
	    Value(make_date(5, January, 2014), 5.15),
	};
	int failures = 0;
	std::unique_ptr<TimeSeries> ts = std::make_unique<TimeSeries>(std::end(values) - std::begin(values), values);
	ts->add(std::end(additional_values) - std::begin(additional_values), additional_values);
	IndexId id = make_index_id(IsdaIndex::CAD_BA_CDOR, TENOR_12M);
	FixingDataService fixings;
	fixings.set_fixings(id, std::move(ts));
	TimeSeries *ts1 = fixings.get_fixings(id);
	if (ts1 == nullptr)
		return 1;
	for (auto v = std::begin(values); v != std::end(values); v++) {
		double result = 0;
		if (!ts1->find(v->date(), result))
			failures++;
		else
			failures += (v->value() != result);
	}
	for (auto v = std::begin(additional_values); v != std::end(additional_values); v++) {
		double result = 0;
		if (!ts1->find(v->date(), result))
			failures++;
		else
			failures += (v->value() != result);
	}
	double result = 0;
	if (ts1->find(make_date(28, December, 2013), result))
		failures++;
	if (ts1->find(make_date(4, January, 2014), result))
		failures++;
	if (ts1->find(make_date(6, January, 2014), result))
		failures++;
	return failures;
}

int test_fixings()
{
	int failure_count = 0;
	failure_count += test_timeseries();
	if (failure_count != 0)
		printf("Fixing Tests FAILED\n");
	else
		printf("Fixing Tests OK\n");
	return failure_count;
}

} // namespace redukti
