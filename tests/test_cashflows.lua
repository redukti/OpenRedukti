-- DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
--
-- Contributor(s):
--
-- The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
-- The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
-- Authors: Dibyendu Majumdar
--
-- Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
--
-- The contents of this file are subject to the the GNU General Public License
-- Version 3 (https://www.gnu.org/licenses/gpl.txt).


local date, date_parts = redukti.date, redukti.date_parts
local index = redukti.index
local utils = assert(require('utils'))
local deposit_rate = utils.deposit_rate

function test_cashflow_build()

	local t = 
	{ -- Collection
		-- stream
		{
			-- cashflow
			{
				type = 'simple',
				currency = 'USD',
				amount = 40523611.1111111,
				payment_date = date(6, 7, 2017)
			},
			-- cashflow
			{
				type = 'ois',
				accrual_start_date = date(5, 7, 2016),
				accrual_end_date = date(3, 7, 2017),
				notional = 815000000,
				index = 'USD-Federal Funds-H.15-OIS-COMPOUND',
				day_count_fraction = 'ACT/360',
				payment_date = date(6, 7, 2017)
			},
			-- cashflow
			{
				type = 'floating',
				compounding_method = 'FLAT',
				currency = 'USD',
				day_count_fraction = 'ACT/360',
				payment_date = date(26,10,2017),
				periods = {
					{
						notional = 10000000,
						accrual_start_date = date(26, 7, 2017),
						accrual_end_date = date(28, 8, 2017),
						index = 'USD-LIBOR-BBA',
						tenor = '1M',
						spread = -0.0009
					},
					{
						notional = 10000000,
						accrual_start_date = date(28, 8, 2017),
						accrual_end_date = date(26, 9, 2017),
						index = 'USD-LIBOR-BBA',
						tenor = '1M',
						spread = -0.0009
					},
					{
						notional = 10000000,
						accrual_start_date = date(26, 9, 2017),
						accrual_end_date = date(26, 10, 2017),
						index = 'USD-LIBOR-BBA',
						tenor = '1M',
						spread = -0.0009
					}
				}
			}
		},
		-- stream
		{
			-- cashflow
			{
				type = 'fra',
				currency = 'USD',
				day_count_fraction = 'ACT/360',
				payment_date = date(8, 9, 2014),
				notional = 15000000,
				fixed_rate = 0.015,
				accrual_start_date = date(8, 9, 2014),
				accrual_end_date = date(20, 11, 2014),
				index = 'USD-LIBOR-BBA',
				tenor = '2M',
				index2 = 'USD-LIBOR-BBA',
				tenor2 = '3M',
			}
		}
	}

	local flows = redukti.cashflows(t)
	local expected = [[
streams { cashflows { simple { currency: USD amount: 40523611.1111111 payment_date: 42922 } } cashflows { ois { index: USD_Federal_Funds_H_15_OIS_COMPOUND notional: 815000000 accrual_start_date: 42556 accrual_end_date: 42919 payment_date: 42922 day_count_fraction: ACT_360 } } cashflows { floating { currency: USD floating_periods { notional: 10000000 spread: -0.0009 accrual_start_date: 42942 accrual_end_date: 42975 index: USD_LIBOR_BBA tenor: TENOR_1M } floating_periods { notional: 10000000 spread: -0.0009 accrual_start_date: 42975 accrual_end_date: 43004 index: USD_LIBOR_BBA tenor: TENOR_1M } floating_periods { notional: 10000000 spread: -0.0009 accrual_start_date: 43004 accrual_end_date: 43034 index: USD_LIBOR_BBA tenor: TENOR_1M } compounding_method: FLAT day_count_fraction: ACT_360 payment_date: 43034 } } } streams { cashflows { fra { currency: USD fixed_rate: 0.015 payment_date: 41890 floating_period { notional: 15000000 accrual_start_date: 41890 accrual_end_date: 41963 index: USD_LIBOR_BBA tenor: TENOR_2M index2: USD_LIBOR_BBA tenor2: TENOR_3M } day_count_fraction: ACT_360 } } }]]
	local got = tostring(flows)
	--print(got)
	assert(got == expected)
end


function test_deposit_rate()
	local flows = deposit_rate(date(25, 7, 2013), 'USD', 'FEDFUND', '1Y', 0.00140)
	local expected = [[
streams { cashflows { simple { currency: USD amount: -1000000 payment_date: 41484 } } cashflows { simple { currency: USD amount: 1001419.4444444444 payment_date: 41849 } } }]]
	local got = tostring(flows)
	assert(got == expected)
end

test_cashflow_build()
test_deposit_rate()

