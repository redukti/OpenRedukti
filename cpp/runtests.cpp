/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti.
 * The Initial Developer of the Original Software is REDUKTI LIMITED.
 *
 * Portions Copyright 2016-2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 2 (http://www.gnu.org/licenses/old-licenses/gpl-2.0.html).
 */
#include <allocators.h>
#include <autodiff.h>
#include <calendars.h>
#include <datasource.h>
#include <date.h>
#include <dayfractions.h>
#include <hashtable.h>
#include <index.h>
#include <statistics.h>
#include <converters.h>
#include <fixings.h>
#include <schedule.h>
#include <interpolators.h>
#include <cashflow_pricing.h>
#include <bootstrap.h>

using namespace redukti;

int main() {
	int rc = 0;
	rc += test_allocators();
	rc += test_date();
	rc += test_hashtable();
	rc += test_calendars();
	rc += test_dayfractions();
	rc += test_index();
	rc += test_autodiff();
	rc += test_datasource();
	rc += test_statistics();
	rc += test_conversions();
	rc += test_fixings();
	rc += test_schedule_generation();
	rc += test_interpolators();
	rc += test_pricing();
	rc += test_bootstrap();
	if (rc != 0)
		printf("FAILED\n");
	else
		printf("OK\n");
	return rc != 0 ? 1 : 0;
}
