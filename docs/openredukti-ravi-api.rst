===============================
OpenRedukti Scripting With Ravi
===============================

OpenRedukti comes with a scripting language `Ravi <http://ravilang.org>`_ that is a dialect of `Lua <http://www.lua.org>`_.
It is much easier to play with and understand how OpenRedukti works by running small examples in the scripting language.

Once you have built OpenRedukti, you will find a commandline program called ``dylan``. To start the scripting environment 
just invoke this command. It works the same way as the Lua command line.

To exit the scripting environment enter ^Z on Windows, and ^D on Unix systems.

The scripting API is a subset of the OpenRedukti C++ API, therefore not every feature is available.

Date Types
==========

Dates are represented as integer values. The api for dates is as follows:

Converting to Date type
-----------------------
reduki.date(day,month,year) 
	takes day, month, year and returns a date
reduki.date(str) 
	Parse a string representation of date. It will detect seperator character '/' or '-'. The formats acceptable are yyyy/mm/dd, dd/mm/yyyy, yyyy-mm-dd, or dd-mm-yyyy
reduki.date{ d,m,y }
	takes table and converts to date
reduki.date() 
	returns today's date

Examples::

	local t = {}
	t.effective_date = redukti.date{1,1,2016}
	t.termination_date = redukti.date{1,1,2017}


Decomposing a Date type
-----------------------
redukti.date_parts(date)
	takes a date and returns day, month, year

Examples::

	function print_date(d)
	    local d,m,y = redukti.date_parts( d )
	    return string.format("%02d/%02d/%d", d, m, y)
	end

	print_date(redukti.date('2017/10/04'))

Date arithmetic
---------------
As a date value represents the number of days since 1899/12/31, you can simply add or subtract from it to get to another date.
Following api is provided for adding a period to date.

redukti.addperiod(date, period)
	Adds period to date and returns resulting date
	Periods are enitities like 1D, 4W, 1M or 15Y.

Examples::

	date = redukti.date()
	print(redukti.addperiod(date, '5Y'))


Holiday Calendars
=================

The api for obtaining and working with calendars is as follows:

redukti.calendar(businesscenters)
	Returns a calendar object. Suported business centers are AUSY, USNY, GBLO, EUTA, JPTO and BRSP. You can specify more than one 
	center if separate the values by a comma.

calendar:advance(date, period [, business_day_convention])
	Advances the given date by the given period and adjusts it to be on a business day. Business day convention defaults to FOLLOWING. 

calendar:advance(date, days [, business_day_convention])
	Moves the date forward by specified number of business days. Business day convention defaults to FOLLOWING.

Examples::

	cal = redukti.calendar('USNY,GBLO')	-- joint calendar
	today = redukti.date()	-- today's date
	month_ago = cal:advance(today, "-1M", "FOLLOWING")
	month_hence = cal:advance(today, "1M", "FOLLOWING")
	print(today, month_ago, month_hence)

Day Count Fractions
===================

redukti.dayfraction(name)
	Returns a Day Count Fraction object for the specified name - name must be ISDA / FpML compatible.

dayfraction:fraction(d1, d2)
	Calculates the day fraction between two dates.

Examples::

	fraction = redukti.dayfraction('ACT/360')
	today = redukti.date()
	two_years_later = redukti.addperiod(today, "2Y")
	print(fraction:fraction(today, two_years_later))


Indices
=======

redukti.index(currency, index_family, tenor)
	Returns an index object

redukti.index(isda_index, tenor)
	Returns an index object

index:maturity(value_date)
	Returns maturity date given the value date, as per the index

index:adjust_date(date, days)
	Adds/subtracts the number of business days specified and adjusts the resulting date as per fixing calendar

index:date_components(date)
	Computes fixing date, value date and maturity date given the start date and returns all three

Examples::

	idx = redukti.index('USD', 'LIBOR', '1W')
	dt = redukti.date(23,10,2016)
	adjusted = idx:adjust_date(dt, 1)
	assert(adjusted == redukti.date(24,10,2016))
	fixing_dt, value_dt, maturity_dt = idx:date_components(adjusted)
	assert(maturity_dt == redukti.date(31,10,2016))
	assert(value_dt == adjusted)

Automatic Differentiation
=========================

redukti.adouble1{ n1, n2, ... }
	Returns an adouble object for each array value. The assumption is that the array values are part of a multivariate function, and therefore each value is treated as variable. The total number of variables is equal to the numeber of values in the array. The adouble objects are set up for first order derrivative computation. Note that the maximum number of allowed variables is 100 to keep memory usage in check. 

redukti.adouble1(n1, n2, ...)
	Same as above but input is in the form of parameters rather than an array.

redukti.adouble2{ n1, n2, ... }
	As above, but the returned adouble objects are set up to compute second order derivatives too.

redukti.adouble2(n1, n2, ...)
	Same as above but input is in the form of parameters rather than an array.

adouble:gradient()
	Returns the gradient as an array

adouble:hessian()
	Only available if second order derivatives are being calculated. The hessian is returned as a table of arrays where each row is an array.

adouble:abs()
	Returns the absolute value

abouble:pow(n)
	Returns adouble object raised to power n

adouble:sqrt()
	Return square root of adouble object

adouble:exp()
	Returns exp(adouble)

adouble:log()
	Returns log(adouble), where log is natural logarithm

adouble:cos()
	Returns cos(adouble)

adouble:sin()
	Return sin(adouble)

adouble:tan()
	Returns tan(adouble)

Examples::

	-- compute derivate of x^2, where x = 5

	x = redukti.adouble1(5.0)
	ans = x:pow(2)
	print(ans)

Above outputs::

	{
	  value=25.0,
	  firstorder = {
	    [1] = 10.0
	  }
	}

This tells you that the value of x^2 is 25.0, and derivative is 10.0 - i.e. 2*x, as you would expect.

Here is another example::

	x, y, z = redukti.adouble2 {5.0, 3.0, 6.0}
	-- compute x + y + z
	added = x + y + z
	multiplied = x * y * z

	print(added)
	print(multiplied)

Results in following output::

	{
	  value=14.0,
	  firstorder = {
	    [1] = 1.0,
	    [2] = 1.0,
	    [3] = 1.0
	  }
	  secondorder = {
	  }
	}
	{
	  value=90.0,
	  firstorder = {
	    [1] = 18.0,
	    [2] = 30.0,
	    [3] = 15.0
	  }
	  secondorder = {
	    [1, 2] = 6.0,
	    [1, 3] = 3.0,
	    [2, 1] = 6.0,
	    [2, 3] = 5.0,
	    [3, 1] = 3.0,
	    [3, 2] = 5.0
	  }
	}

Since the matrix data for a second order derivative can grow very large, the scripting api restricts the number of allowed variables in a single object to 100.

Another example::

	x, y = redukti.adouble2 {6, 7}
	ans = x:pow(2) * y:pow(2)
	print(ans)

Results in::

	{
	  value=1764.0,
	  firstorder = {
	    [1] = 588.0,
	    [2] = 504.0
	  }
	  secondorder = {
	    [1, 1] = 98.0,
	    [1, 2] = 168.0,
	    [2, 1] = 168.0,
	    [2, 2] = 72.0
	  }
	}

Calculation Schedules
=====================

redukti.schedule { parameters }
	Builds a calculation schedule and returns 3 arrays : adjusted start dates, adjusted end dates, adjusted payment dates. Note that some payment dates may be set to 0 - this means that there is no payment in that calculation period. For instance when compounding, or in zero coupon streams, payments do not occur with every period. 

	The parameters may include following:

	effective_date
		(required) unadjusted effective date - this defines the start of the schedule

	termination_date
		(required if term not present) unadjusted termination date - this defines the end date of the schedule.

	term
		(required if termination_date not present) term is the length of the transaction, e.g. 5Y.

	payment_frequency
		(required) the frequency of payment, may be 1T for single payment, or supported tenor values upto 1Y.

	payment_day_convention
		BusinessDayConvention to be used to adjust payment dates

	payment_business_centers
		Business Centers to be use for computing payment holidays

	roll_convention
		FpML defined RollConvention for deciding how to calculate period start/end dates

	calculation_frequency
		the frequency of calculating accruals, also equal to the frequency at which resets occur for non-OIS streams. Must be <= payment frequency

	calculation_day_convention
		the BusinessDayConvention to be used for adjusting calculation period dates

	calculation_business_centers
		Business Centers to be used for computing holidays when adjusting calculation period dates

	first_regular_period_start_date
		unadjusted start date of the first regular period, implies front stub

	last_regular_period_end_date
		unadjusted end date of the last regular period, implies back stub

	first_payment_date
		used if first payment does not occur at the first possible payment date

	last_regular_payment_date 
		used if last regular period payment date does not occur at the last possible regular payment date

Example::

	t = {}
	t.effective_date = redukti.date{1,1,2016}
	t.termination_date = redukti.date{1,1,2017}
	t.payment_frequency = "3M"
	t.payment_business_centers = "GBLO,USNY"
	t.payment_day_convention = "MODFOLLOWING"

	starts, ends, pays = redukti.schedule(t)

	function print_date(d)
	    local d,m,y = redukti.date_parts( d )
	    return string.format("%02d/%02d/%d", d, m, y)
	end

	for i=1,#starts do
		print(i, '  adjusted start date ' .. print_date(starts[i]))
		print(i, '    adjusted end date ' .. print_date(ends[i]))
		if pays[i] ~= 0 then
			print(i, 'adjusted payment date ' .. print_date(ends[i]))
		else
			print(i, '           no payment ')
		end			
	end

This outputs::

	1         adjusted start date 04/01/2016
	1           adjusted end date 01/04/2016
	1       adjusted payment date 01/04/2016
	2         adjusted start date 01/04/2016
	2           adjusted end date 01/07/2016
	2       adjusted payment date 01/07/2016
	3         adjusted start date 01/07/2016
	3           adjusted end date 03/10/2016
	3       adjusted payment date 03/10/2016
	4         adjusted start date 03/10/2016
	4           adjusted end date 03/01/2017
	4       adjusted payment date 03/01/2017

Another example::

	t = {}
	t.effective_date = redukti.date{25,11,2016}
	t.termination_date = redukti.date{25,11,2017}
	t.payment_frequency = "3M"
	t.payment_business_centers = "GBLO,USNY"
	t.payment_day_convention = "MODFOLLOWING"
	t.calculation_frequency = "1M"
	t.calculation_business_centers = "GBLO,USNY"
	t.calculation_day_convention = "MODFOLLOWING"

	starts, ends, pays = redukti.schedule(t)

This time the output is::

	1         adjusted start date 25/11/2016
	1           adjusted end date 28/12/2016
	1                  no payment
	2         adjusted start date 28/12/2016
	2           adjusted end date 25/01/2017
	2                  no payment
	3         adjusted start date 25/01/2017
	3           adjusted end date 27/02/2017
	3       adjusted payment date 27/02/2017
	4         adjusted start date 27/02/2017
	4           adjusted end date 27/03/2017
	4                  no payment
	5         adjusted start date 27/03/2017
	5           adjusted end date 25/04/2017
	5                  no payment
	6         adjusted start date 25/04/2017
	6           adjusted end date 25/05/2017
	6       adjusted payment date 25/05/2017
	7         adjusted start date 25/05/2017
	7           adjusted end date 26/06/2017
	7                  no payment
	8         adjusted start date 26/06/2017
	8           adjusted end date 25/07/2017
	8                  no payment
	9         adjusted start date 25/07/2017
	9           adjusted end date 25/08/2017
	9       adjusted payment date 25/08/2017
	10        adjusted start date 25/08/2017
	10          adjusted end date 25/09/2017
	10                 no payment
	11        adjusted start date 25/09/2017
	11          adjusted end date 25/10/2017
	11                 no payment
	12        adjusted start date 25/10/2017
	12          adjusted end date 27/11/2017
	12      adjusted payment date 27/11/2017

Note that in the examples above I did not specify roll convention so this was inferred. 

Interpolators
=============

redukti.interpolator { parameters }
	Returns an interpolator object

	Parameters can be following:

	x
		values to be used for x axis

	y
		values to be used for y axis (interpolation occurs in y axis)

	order
		1 if first order derivatives are needed; 2 will generate second order derivatives also

	interpolator
		The type of interpolator, e.g. Linear, LogLinear, CubicSplineNatural, LogCubicSplineNatural, MonotoneConvex

interpolator:interpolate(x [, n])
	Returns the interpolated value. If the optional parameter n is 1 the return value is an adouble containing the interpolated value as well as the derivatives computed using automatic differentiation. If the optional parameter n is 2 then the return value is an adouble but computed using numeric differentiation. The latter is for testing purposes only.

Examples::

	x = {0.01, 0.02, 0.03, 0.04, 0.05}

	y = {1000000.0, 20004.0, 300000.5, 4000000.0, 900000.0}

	interp1 = redukti.interpolator {
		x = x,
		y = y,
		interpolator = 'CubicSplineNatural',
		order = 2
	}

	interp2 = redukti.interpolator {
		x = x,
		y = y,
		interpolator = 'MonotoneConvex',
		order = 2
	}

	print(interp1:interpolate(0.035, 1))

	print(interp2:interpolate(0.035, 1))

Interest Rate Curves
====================

redukti.curve { parameters }
	Sets up an interest rate curve. The parameters are as follows.

	reference_date
		The date of the curve, all maturities are with respect to this date

	maturities
		An array of maturity dates

	values
		An array of numbers representing zero rates or discount factors

	value_type 
		'ZeroRate' or 'DiscountFactor' - indicates the type of value

	interpolator
		Interpolator type. For ZeroRate curves, use Linear, CubicSplineNatural, MonotoneConvex. For discount factor curves, use LogLinear, LogCubicSplineNatural

	currency
		Currency, forms part of curve's id

	index_family
		IndexFamily, forms part of curve's id

	tenor
		Curve's tenor, forms part of curve's id

	order
		If 1 first order derivatives will be computed, if 2 additionally second order derivatives will be computed

	curve_type
		'Forward' or 'Discount' - this is a logical marker for how the curve will be used, forms part of curve's id 

curve:zero_rate(date)
	Returns the zero rate 

curve:discount_factor(date)
	Returns the discount factor

curve:sensitivities(date)
	Computes and returns sensitivities for given date with respect to the curve pillars.

curve:values()
	Returns three arrays - maturies, zero rates and discount factors

Examples:

As creating a curve manually is tedious, often it is easier to import data from files. For an example of this please see the the example in `test_zerocurve.lua <https://github.com/redukti/OpenRedukti/blob/master/tests/test_zerocurve.lua>`_.

Time Series / Fixings
=====================

redukti.fixing_service { data }
	Creates a FixingService object. The data must be a table containing values indexed by Index type. Each value must be a table that has fixings indexed by date. That is, of the form:

::

	{
		[index1] = {
			[date1] = index1value1,
			[date2] = index1value2
		},
		[index2] = {
			[date1] = index2value1,
			[date2] = index2value2
		}
	}


fixing_service:fixing(index, date)
	Returns the fixing for given index and date. If not found then returns nil

Examples:

This too is easier to load from files, hence I will refer to the following script as an example.

Cashflows
=========

redukti.cashflows { data }
	Sets up a CFCollection object and returns it.

The contents of ``{ data }`` mirrors the CFCollection protocol buffers type, except that it is expressed as a Lua table. Here is an example::

	t = 
	{ -- Collection
		-- first stream
		{
			-- cashflow
			{
				type = 'simple',
				currency = 'USD',
				amount = 40523611.1111111,
				payment_date = redukti.date(6, 7, 2017)
			},
			-- cashflow
			{
				type = 'ois',
				accrual_start_date = redukti.date(5, 7, 2016),
				accrual_end_date = redukti.date(3, 7, 2017),
				notional = 815000000,
				index = 'USD-Federal Funds-H.15-OIS-COMPOUND',
				day_count_fraction = 'ACT/360',
				payment_date = redukti.date(6, 7, 2017)
			},
			-- cashflow
			{
				type = 'floating',
				compounding_method = 'FLAT',
				currency = 'USD',
				day_count_fraction = 'ACT/360',
				payment_date = redukti.date(26,10,2017),
				periods = {
					{
						notional = 10000000,
						accrual_start_date = redukti.date(26, 7, 2017),
						accrual_end_date = redukti.date(28, 8, 2017),
						index = 'USD-LIBOR-BBA',
						tenor = '1M',
						spread = -0.0009
					},
					{
						notional = 10000000,
						accrual_start_date = redukti.date(28, 8, 2017),
						accrual_end_date = redukti.date(26, 9, 2017),
						index = 'USD-LIBOR-BBA',
						tenor = '1M',
						spread = -0.0009
					},
					{
						notional = 10000000,
						accrual_start_date = redukti.date(26, 9, 2017),
						accrual_end_date = redukti.date(26, 10, 2017),
						index = 'USD-LIBOR-BBA',
						tenor = '1M',
						spread = -0.0009
					}
				}
			}
		},
		-- second stream
		{
			-- cashflow
			{
				type = 'fra',
				currency = 'USD',
				day_count_fraction = 'ACT/360',
				payment_date = redukti.date(8, 9, 2014),
				notional = 15000000,
				fixed_rate = 0.015,
				accrual_start_date = redukti.date(8, 9, 2014),
				accrual_end_date = redukti.date(20, 11, 2014),
				index = 'USD-LIBOR-BBA',
				tenor = '2M',
				index2 = 'USD-LIBOR-BBA',
				tenor2 = '3M',
			}
		}
	}

	flows = redukti.cashflows(t)
	print(tostring(flows))

When you run this the output will dump the generated protocol buffers value in JSON like format.

Utility for Loading Data
========================

The scripting api contains following utility for loading data from CSV files.

redukti.loadcsv { parameters }
	Loads data from a CSV file. Returns a table where each element represents a row in the CSV file.

	Following parameters are supported:

	file
		Specifies the path and name of the file to read from

	conversion
		Contains a string where each character represents a conversion rule for a column. Following rules are allowed:

		's'
			Intepret as string field

		'n'
			Interpret as number field

		'i'
			Interpret as integer field

		'd'
			Interpret as date field

		'-'
			Ignore column, field set to nil

	headings
		A value of true means that the source file has headings

	fields
		A value of true means that each field will be named by the column heading

Suppose that a file contains::

	index,tenor,date,fixing
	EUR-EONIA-OIS-COMPOUND,1D,02/01/2007,0.036
	EUR-EURIBOR-Reuters,1W,02/01/2007,0.03614
	EUR-EURIBOR-Reuters,2W,02/01/2007,0.03615
	EUR-EURIBOR-Reuters,1M,02/01/2007,0.03629
	EUR-EURIBOR-Reuters,3M,02/01/2007,0.03725
	EUR-EURIBOR-Reuters,6M,02/01/2007,0.03857
	EUR-EURIBOR-Reuters,12M,02/01/2007,0.0403
	EUR-EONIA-OIS-COMPOUND,1D,03/01/2007,0.036

Then we can load this using following::

	fixings = redukti.loadcsv { file=filename, conversion='ssdn', heading=true, fields=true }

Here is a dump of the first three lines of the table::

	> table_print(fixings)
	[1] => table
	    (
	       [tenor] => 1D
	       [index] => EUR-EONIA-OIS-COMPOUND
	       [fixing] => 0.036
	       [date] => 39084
	    )
	[2] => table
	    (
	       [tenor] => 1W
	       [index] => EUR-EURIBOR-Reuters
	       [fixing] => 0.03614
	       [date] => 39084
	    )
	[3] => table
	    (
	       [tenor] => 2W
	       [index] => EUR-EURIBOR-Reuters
	       [fixing] => 0.03615
	       [date] => 39084
	    )

Building Curves
===============

redukti.build_curves( business_date, curve_definitions, par_rates [, pricing_script] )
	Builds a set of Zero Rate curves from par rates. Requires curve definitions and par rates as input.

	curve_definitions
		The curve definitions must be presented as a table that mirrors the IRCurveDefinition protocol buffers message type. Values must be keyed by the IRCurveDefinition.id. Each value must be a table conatining the fields corresponding to an IRCurveDefinition. These are:

		group
			CurveGroup

		curve_id
			ID of the curve

		interpolator
			Interpolator to be used, e.g. Linear

		currency
			Currency of the curve

		index_family
			IndexFamily

		interpolated_on
			'ZeroRate' or 'DiscountFactor'

		curve_tenor
			Tenor of the curve

		maturity_generation_rule
			Specifies how maturities for instruments are derived

		tenors
			Optional list of tenor values used if maturities are derived from fixed tenors

	par_rates
		The input par rates must be presented as a table. The elements in the table must have following fields.

		curve_id
			a IRCurveDefinition.id that is defined in the curve definitions

		instrument_type
			The type of instrument. This is mapped to a Ravi / Lua scripting function name.

		instrument_id
			The id of the instrument - this has to be in a specific format.

		par_rate
			The par rate

		forward_curve_id
			IRCurveDefinition.id of the curve to be used for forward rates

		discount_curve_id
			IRCurveDefinition.id of the curve to be used for discounting

		floating_tenor
			The tenor to be used on floating leg of the instrument.

	pricing_script
		This is a Lua script that is responsible for generating the cashflow structure for each instrument used in the curve. See the default script named ``pricing.lua`` that is supplied with OpenRedukti. If a name isn't supplied the script defaults to the one used when building OpenRedukti; this is useful for testing but not very good for deployment as the path to the script is baked in at compile time! 

Examples:

Please see the function ``build_curves`` in `utils.lua <https://github.com/redukti/OpenRedukti/blob/master/tests/utils.lua>`_.


Cashflow Pricing
================

This is complex process involving several steps.

1. First of all you need a set of zero curves, either bootstrapped from par rates, or obtained from another source.
2. You need to setup a curve mapper.
3. You need fixings.
4. You need to setup a ValuationContext.
5. You then convert the Cashflows to an internal format ready for pricing.
6. Next you need to setup a curve provider.
7. Finally you invoke the pricing function to compute the PV and the sensitivities of the cashflows.

Curve Mapper
------------

The purpose of the Curve Mapper is to allow mapping of logical curves, so that the cashflow generator can reference curves without knowing how these curves will actually be delivered. So this provides a level of indirection.

The api is described below.

redukti.curve_mapper()
	Sets up a curve mapper object

redukti.pricing_curve_id { parameters }
	Returns a logical curve id. The parameters are:

	curve_type
		'Forward' or 'Discount'

	currency
		Currency

	index_family
		IndexFamily

	tenor
		Optional tenor of the curve

curve_mapper:add_mapping( from_id, to_id )
	This sets up a mapping from 'from_id' to 'to_id'. Both ids must be pricing_curve_ids.

Example::

	curve_mapper = redukti.curve_mapper()
	f_eonia_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EONIA'}
	d_eonia_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR' }
	d_euribor_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR'}
	f_euribor_1m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '1M'}
	f_euribor_3m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '3M'}
	f_euribor_6m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '6M'}
	f_euribor_12m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '12M'}
	d_euribor_1m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '1M'}
	d_euribor_3m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '3M'}
	d_euribor_6m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '6M'}
	d_euribor_12m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '12M'}
	f_euribor_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR'}

	-- map request for forward EONIA curve to discount curve
	curve_mapper:add_mapping( f_eonia_id, d_eonia_id )
	-- map any request for EURIBOR Discount to EONIA discount
	curve_mapper:add_mapping( d_euribor_id, d_eonia_id )
	-- map euribor 1m to generic
	curve_mapper:add_mapping( f_euribor_1m_id, f_euribor_id )
	-- map euribor 12m to generic
	curve_mapper:add_mapping( f_euribor_1m_id, f_euribor_id )

ValuationContext
----------------

redukti.valuation_context( parameters, fixing_service)
	Sets up a ValuationContex object. 

	parameters
		Table containing parameters. Parameters are:
	
		reference_date
			Business date

		order
			If 1 first order derivatives will be computed, if 2, second order derivatives will be computed as well.

	fixing_service
		A FixingService object

Example::

	utils = assert(require('utils'))
	fixing_service = utils.load_fixings('../testdata/20121211/fixings.csv')
	business_date = redukti.date('2012/12/11')
	valuation_context = redukti.valuation_context({ reference_date = business_date, order = 1 }, fixing_service)

Cashflow Conversion
-------------------
redukti.prepare_cashflows_for_pricing(valuation_context, curve_mapper, cashflows)
	Converts the supplied cashflows to a format suitable for pricing.

Example::

	flows = deposit_rate(redukti.date(25, 7, 2013), 'EUR', 'EURIBOR', '1Y', 0.00140)
	pricing_cashflows = redukti.prepare_cashflows_for_pricing(valuation_context, curve_mapper, flows)


CurveProvider
-------------
The CurveProvider resolves logical curve ids and maps these to actual curve objects.

redukti.curve_provider()
	Creates a curve_provider object.

curve_provider:add_mapping(pricing_curve_id, curve)
	Adds a mapping from a pricing_curve_id to a Zero Curve object.

Example:

Please see the script `test_pricing.lua <https://github.com/redukti/OpenRedukti/blob/master/tests/test_pricing.lua>`_.

Calculate NPV and Sensitivities
-------------------------------
Finally you can invoke:

redukti.present_value(valuation_context, pricing_cashflows, curve_provider)
	Computes the NPV and sensitivities and returns a PricingResult object.

pricing_result:ok()
	Tells if you pricing was successful

pricing_result:curve_ids()
	Returns an array of curve identifiers used for pricing

pricing_result:delta(curve_id)
	Given a curve id, returns the first order sensitivities as an table of values keyed by curve pillars.

Example:

Please see the script `test_pricing.lua <https://github.com/redukti/OpenRedukti/blob/master/tests/test_pricing.lua>`_.
