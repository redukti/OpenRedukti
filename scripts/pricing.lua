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

local bond_templates = {
	USD_GOVT_BOND_6M = {
		payment_calendar = "USNY",
		start_delay = 2,
		fixed_day_fraction = "ACT/ACT.ISMA",
		payment_day_convention = "MODFOLLOWING",
		fixed_payment_frequency = "6M"
	}
}

local deposit_templates = {
  USD_FEDFUND = {
    payment_calendar = "USNY",
  },
  USD_LIBOR = {
    payment_calendar = "USNY",
  },
  USD_LIBOR_3M = {
    payment_calendar = "USNY",
  },
  USD_LIBOR_6M = {
    payment_calendar = "USNY",
  },
  EUR_EONIA = {
    payment_calendar = "EUTA",
  },
  EUR_EURIBOR = {
    payment_calendar = "EUTA",
  },
  EUR_EURIBOR_3M = {
    payment_calendar = "EUTA",
  },
  EUR_EURIBOR_6M = {
    payment_calendar = "EUTA",
  },
  USD_GOVT_BOND = {
    payment_calendar = "USNY",
	start_delay = 0,
	fixed_day_fraction = 'ACT/365.FIXED',
  }
}

local swap_templates = {
  USD_FEDFUND = {
    start_delay = 2,
    payment_calendar = "USNY",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "6M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "USD-Federal Funds-H.15-OIS-COMPOUND"
  },
  USD_LIBOR = {
    start_delay = 2,
    payment_calendar = "USNY",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "USD-LIBOR-BBA"
  },
  USD_LIBOR_3M = {
    start_delay = 2,
    payment_calendar = "USNY",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "USD-LIBOR-BBA"
  },
  USD_LIBOR_6M = {
    start_delay = 2,
    payment_calendar = "USNY",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "6M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "USD-LIBOR-BBA"
  },
  EUR_EONIA = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "1Y",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EONIA-OIS-COMPOUND"
  },
  EUR_EURIBOR = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EURIBOR-Reuters"
  },
  EUR_EURIBOR_3M = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EURIBOR-Reuters"
  },
  EUR_EURIBOR_6M = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "6M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EURIBOR-Reuters"
  }
}

local fra_templates = {
  USD_LIBOR_3M = {
    start_delay = 2,
    payment_calendar = "USNY",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "USD-LIBOR-BBA"
  },
  USD_LIBOR_6M = {
    start_delay = 2,
    payment_calendar = "USNY",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "6M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "USD-LIBOR-BBA"
  },
  EUR_EURIBOR = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EURIBOR-Reuters"
  },
  EUR_EURIBOR_3M = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "3M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EURIBOR-Reuters"
  },
  EUR_EURIBOR_6M = {
    start_delay = 2,
    payment_calendar = "EUTA",
    fixed_payment_frequency = "1Y",
    floating_payment_frequency = "6M",
    fixed_day_fraction = "ACT/360",
    floating_day_fraction = "ACT/360",
    payment_day_convention = "MODFOLLOWING",
	floating_index = "EUR-EURIBOR-Reuters"
  }
}


function deposit_rate(business_date: integer, ccy, index_family, tenor, inst_id, rate: number)
	local template = deposit_templates[ccy .. '_' .. index_family .. '_' .. tenor]
	if not template then
		template = deposit_templates[ccy .. '_' .. index_family]
	end
	if not template then
		error("template for not found for " .. ccy .. '_' .. index_family .. '_' .. tenor)
		return nil
	end
	local day_count = 'ACT/365.FIXED' 
	local start_delay = 0
	if template.start_delay then
		start_delay = template.start_delay
	else
		if (ccy ~= 'GBP') then
			if (inst_id == '1D') then
				start_delay = 0
			elseif (inst_id == '2D') then
				start_delay = 1
			else
				start_delay = 2
			end
		end
	end
	if template.fixed_day_fraction then
		day_count = template.fixed_day_fraction
	else
		if (ccy ~= 'GBP') then
			day_count = 'ACT/360'
		end
	end
	local calendar = redukti.calendar(template.payment_calendar)
	if not calendar then
		error("payment calendar not found for " .. template.payment_calendar)
		return nil
	end

	--local index = redukti.index(ccy, index_family, inst_id)
	
	local notional: number = 1000000.0
	local start_date: integer = calendar:advance(business_date, start_delay)
	--start_date  = index:adjust_date(business_date, start_delay)
	local maturity_date: integer = calendar:advance(start_date, inst_id)
	--maturity_date = index:maturity(start_date)
	local dayfraction = redukti.dayfraction(day_count)
	if not dayfraction then
		return nil
	end
	local tau: number = dayfraction:fraction(start_date, maturity_date)
	local cf = 
	{
		-- stream
		{
			-- cashflow
			{
				type = 'simple',
				currency = ccy,
				amount = -notional,
				payment_date = start_date
			},
			-- cashflow
			{
				type = 'simple',
				currency = ccy,
				amount = (1.0 + rate * tau) * notional,
				payment_date = maturity_date			
			}
		}
	}
	local flows = redukti.cashflows(cf)
    --print(ravitype(flows))
    --print(flows)
	return maturity_date, flows
end

-- price a swap with only one curve
function swap_rate_sc(today: integer, ccy, index_family, tenor, inst_id, par_rate: number)
	-- retrieve a template
	-- if tenor specific template is not found then
	-- look for a generic template for the given index and currency
	local template = swap_templates[ccy .. '_' .. index_family .. '_' .. tenor]
	if not template then
		template = swap_templates[ccy .. '_' .. index_family]
	end
	if not template then
		error("template for not found for " .. ccy .. '_' .. index_family .. '_' .. tenor)
		return nil
	end
  
  	local notional: number = 1000000.0
	local calendar = redukti.calendar(template.payment_calendar)
	if not calendar then
		return nil
	end

    -- swap start date
	local start_date: integer = calendar:advance(today, template.start_delay)
	--local index = redukti.index(ccy, index_family, inst_id)
	--start_date = index:adjust_date(today, template.start_delay)

	-- swap end date
	local maturity_date: integer = redukti.addperiod(start_date, inst_id)

	local S = {}
	S.effective_date = start_date
	S.termination_date = maturity_date
	S.payment_business_centers = template.payment_calendar
	S.payment_day_convention = template.payment_day_convention
	S.payment_frequency = template.fixed_payment_frequency

	local fixstarts: integer[], fixends: integer[], fixpays: integer[] = redukti.schedule(S)
	assert(fixstarts ~= nil)
  
	local daycount = redukti.dayfraction(template.fixed_day_fraction)
	if not daycount then
		return nil
	end

	local fixscalars: number[] = {}
	local i
	for i = 1,#fixstarts do
    	fixscalars[i] = notional * par_rate * 
       		@number daycount:fraction(fixstarts[i], fixends[i])
	end

	local cfstream = {}	
	local cfcollection = { cfstream }

	cfstream[1] = {
		type = 'simple',
		currency = ccy,
		amount = -notional,
		payment_date = start_date
	}
	cfstream[2] = {
		type = 'simple',
		currency = ccy,
		amount = notional,
		payment_date = maturity_date		
	}
	for i = 1,#fixscalars do
		cfstream[i+2] = {
			type = 'simple',
			currency = ccy,
			amount = fixscalars[i],
			payment_date = fixpays[i]		
		}
	end

	local flows = redukti.cashflows(cfcollection)
    --print(flows)
	return maturity_date, flows
end

-- price a vanilla (non compounding) swap with separate forward curve and discount curve
function swap_rate(today: integer, ccy, index_family, tenor, inst_id, par_rate: number, floating_tenor)
	-- retrieve a template
	-- if tenor specifi template is not found then
	-- look for a generic template for the given index and currency
	local template = swap_templates[ccy .. '_' .. index_family .. '_' .. tenor]
	if not template then
		template = swap_templates[ccy .. '_' .. index_family]
	end
	if not template then
		error("template for not found for " .. ccy .. '_' .. index_family .. '_' .. tenor)
		return nil
	end
  
	if not floating_tenor then
		floating_tenor = template.floating_payment_frequency
	end
  	local notional: number = 1000000.0
	local calendar = redukti.calendar(template.payment_calendar)
	if not calendar then
		error("calendar not found for " .. template.payment_calendar)
		return nil
	end

    -- swap start date
	local start_date: integer = calendar:advance(today, template.start_delay)
	-- swap end date
	local maturity_date: integer = redukti.addperiod(start_date, inst_id)

	local S = {}
	S.effective_date = start_date
	S.termination_date = maturity_date
	S.payment_business_centers = template.payment_calendar
	S.payment_day_convention = template.payment_day_convention
	S.payment_frequency = floating_tenor

	local floatstarts: integer[], floatends: integer[], floatpays: integer[] = redukti.schedule(S)

	S.payment_frequency = template.fixed_payment_frequency

	local fixstarts: integer[], fixends: integer[], fixpays: integer[] = redukti.schedule(S)
  
	local fixed_daycount = redukti.dayfraction(template.fixed_day_fraction)
	if not fixed_daycount then
		return nil
	end

	local fixscalars: number[] = {}
	local i
	for i = 1,#fixstarts do
    	fixscalars[i] = notional * par_rate * 
       		@number fixed_daycount:fraction(fixstarts[i], fixends[i])
	end

	-- fixed stream pays and floating stream receives
	local fixedstream = { factor=-1.0 }
	local floatstream = { factor=1.0 }	
	local cfcollection = { fixedstream, floatstream }

	for i = 1,#fixscalars do
		fixedstream[i] = {
			type = 'simple',
			currency = ccy,
			amount = fixscalars[i],
			payment_date = fixpays[i]		
		}
	end
	for i = 1,#floatstarts do
		floatstream[i] = {
			type = 'floating',
			currency = ccy,
			day_count_fraction = template.floating_day_fraction,
			payment_date = floatpays[i],
			periods = {
				{
					notional = notional,
					accrual_start_date = floatstarts[i],
					accrual_end_date = floatends[i],
					index = template.floating_index,
					tenor = floating_tenor,
					spread = 0.0
				}
			}
		}
	end
	local flows = redukti.cashflows(cfcollection)
	return maturity_date, flows
end


function imm_future(today, ccy, index_family, tenor, inst_id, par_rate)
	local notional: number = 1000000.0
	local index = redukti.index(ccy, index_family, "3M")
	local start_date, end_date = redukti.future_imm_dates(ccy, inst_id);

	local fixing_date, value_date, maturity_date = index:date_components(start_date)
	if currency == "AUD" then
		maturity_date = index:date_adjust(start_date, 90)
	elseif currency == "DKK" or currency == "SEK" or currency == "NOK" then
		maturity_date = end_date
	end
  
	local day_count = "ACT/365.FIXED"
	if (currency ~= "GBP") then
		day_count = "ACT/360"
	end  
  
	local dayfraction = redukti.dayfraction(day_count)
	if not dayfraction then
		return nil
	end
  
	local frac_period = dayfraction:fraction(start_date, end_date)
	local frac_index = dayfraction:fraction(value_date, maturity_date)
	local ratio = frac_period / frac_index
	local amount = notional * (1 + par_rate * frac_index)^ratio

	local cf = 
	{
		-- stream
		{
			-- cashflow
			{
				type = 'simple',
				currency = ccy,
				amount = -notional,
				payment_date = start_date
			},
			-- cashflow
			{
				type = 'simple',
				currency = ccy,
				amount = amount,
				payment_date = end_date
			}
		}
	}
	local flows = redukti.cashflows(cf)
	return maturity_date, flows
end

function fra(today, ccy, index_family, tenor, inst_id, par_rate)
	local start_tenor, end_tenor = redukti.fra_tenors(inst_id)
	local fra_term = (end_tenor - start_tenor) .. "M"
	local template = fra_templates[ccy .. '_' .. index_family .. '_' .. fra_term]
	if not template then
		error("template for not found for " .. ccy .. '_' .. index_family .. '_' .. tenor)
		return nil
	end

	local calendar = redukti.calendar(template.payment_calendar)
	local start_date = calendar:advance(today, template.start_delay)
	local fra_start = calendar:advance(start_date, start_tenor .. "M")
	local fra_end = calendar:advance(start_date, end_tenor .. "M")
	local notional = 1000000

	local cf = {
		{
			{
				type = 'fra',
				currency = ccy,
				day_count_fraction = template.floating_day_fraction,
				payment_date = fra_start,
				notional = notional,
				fixed_rate = par_rate,
				accrual_start_date = fra_start,
				accrual_end_date = fra_end,
				index = template.floating_index,
				tenor = fra_term
			}
		}
	}
	local flows = redukti.cashflows(cf)
	return fra_end, flows
end


function test_swap_rate_sc()
	local maturity_date, flows = swap_rate_sc(date(25, 7, 2013), 'USD', 'FEDFUND', '', '4Y', 0.009107)
	local expected_cashflows = [[
streams { cashflows { simple { currency: USD amount: -1000000 payment_date: 41484 } } cashflows { simple { currency: USD amount: 1000000 payment_date: 42945 } } cashflows { simple { currency: USD amount: 9233.4861111111113 payment_date: 41849 } } cashflows { simple { currency: USD amount: 9233.4861111111113 payment_date: 42214 } } cashflows { simple { currency: USD amount: 9258.7833333333328 payment_date: 42580 } } cashflows { simple { currency: USD amount: 9284.0805555555544 payment_date: 42947 } } }]]
	assert(tostring(flows) == expected_cashflows)
end

function test_swap_rate()
	local maturity_date, flows = swap_rate(date(25, 7, 2013), 'USD', 'LIBOR', '3M', '2Y', 0.004888)
	assert(maturity_date == date(29, 7, 2015)) 
	--assert(tostring(flows) == expected_cashflows)
	print(flows)
end

function test_imm_future()
	local maturity_date, flows = imm_future(date(25, 7, 2013), 'USD', 'LIBOR', '3M', 'Dec13', 0.002844)
	--assert(maturity_date == date(29, 7, 2015)) 
	--assert(tostring(flows) == expected_cashflows)
	print(flows)
end

function test_fra()
	local maturity_date, flows = fra(date(25, 7, 2013), 'USD', 'LIBOR', '6M', '12X18F', 0.00645)

	print(flows)
end

function test_deposit_rate()
	local maturity_date, flows = deposit_rate(date(25, 7, 2013), 'USD', 'LIBOR', '3M', '3M', 0.002638)

	print(flows)
end

-- price a bond with par value of 100, and tenor 6M
-- instrument_id must be made up of issue_date:maturity_date:clean_price
-- par_rate must the coupon rate
-- TODO make par value, frequency a parameter
-- TODO and refactor this to a wrapper function
function fixed_bond_100_6M(today: integer, ccy: string, index_family: string, tenor, inst_id: string, coupon_rate: number)

	local issue_date, maturity_date, clean_price_str = string.match(inst_id, "([%d-]+):([%d-]+):([%d.]+)")
	if not issue_date or
		not maturity_date or
		not clean_price_str then
		error("Instrument id is not in correct format")
		return nil
	end

	-- retrieve a template
	-- if tenor specifi template is not found then
	-- look for a generic template for the given index and currency

	local tenor = '6M'
	local template = bond_templates[ccy .. '_' .. index_family .. '_' .. tenor]
	if not template then
		error("template for not found for " .. ccy .. '_' .. index_family .. '_' .. tenor)
		return nil
	end
  
	local units = 10000 -- to make the notional = 1m
	local notional: number = 100.0
	local clean_price: number = tonumber(clean_price_str)
	local fixed_rate: number = coupon_rate
	local calendar = redukti.calendar(template.payment_calendar)
	if not calendar then
		return nil
	end

	-- bond effective date - may be in the past
	local start_date: integer = redukti.date(issue_date)

	-- bond maturity date
	local maturity_date: integer = redukti.date(maturity_date)

	local roll_day, roll_month, roll_year = redukti.date_parts(maturity_date)

	local S = {}
	S.effective_date = start_date
	S.termination_date = maturity_date
	S.payment_business_centers = template.payment_calendar
	S.payment_day_convention = template.payment_day_convention
	S.payment_frequency = template.fixed_payment_frequency
	S.calculation_day_convention = 'N'
	S.stub_type = 'f'
	--S.roll_convention = roll_day >= 30 and 'EOM' or roll_day

	local cutoff_date: integer = calendar:advance(today, template.start_delay)
	local fixstarts: integer[], fixends: integer[], fixpays: integer[], stubtype: string = redukti.schedule(S)
	assert(fixstarts ~= nil)
  
	--for i = 1,#fixstarts do
	-- print("Start " .. fixstarts[i] .. " end " .. fixends[i] .. " pay on " .. fixpays[i] .. ' Cutoff? ' .. (fixpays[i] > cutoff_date and 'N' or 'Y'))
	--end

	local daycount = redukti.dayfraction(template.fixed_day_fraction)
	if not daycount then
		return nil
	end

	local fixscalars: number[] = {}
	local i
	for i = 1,#fixstarts do
		if (template.fixed_day_fraction == 'ACT/ACT.ISMA') then
			local refstart: integer
			local refend: integer
			if i < #fixstarts then
				refstart = redukti.addperiod(fixends[i], '-' .. template.fixed_payment_frequency)
				refend = fixends[i]
			else
				refstart = fixstarts[i]
				refend = redukti.addperiod(fixstarts[i], template.fixed_payment_frequency)
			end
	    	fixscalars[i] = notional * units * fixed_rate * 
    	   		@number daycount:fraction(fixstarts[i], fixends[i], refstart, refend)
			--print("Start " .. fixstarts[i] .. " end " .. fixends[i] .. " refstart " .. refstart .. ' refend ' .. refend)
		else
	    	fixscalars[i] = notional * units * fixed_rate * 
    	   		@number daycount:fraction(fixstarts[i], fixends[i])
    	end
	end

	local cfstream = {}	
	local cfcollection = { cfstream }

	cfstream[1] = {
		type = 'simple',
		currency = ccy,
		amount = -clean_price * units,
		payment_date = today
	}
	cfstream[2] = {
		type = 'simple',
		currency = ccy,
		amount = notional * units,
		payment_date = maturity_date		
	}
	local j = 3
	for i = 1,#fixscalars do
		if fixpays[i] >= cutoff_date then
			cfstream[j] = {
				type = 'simple',
				currency = ccy,
				amount = fixscalars[i],
				payment_date = fixpays[i]		
			}
			j = j + 1
		end
	end

	local flows = redukti.cashflows(cfcollection)
    -- print(flows)
	return maturity_date, flows
end

--test_swap_rate_sc()
--test_swap_rate()
--test_imm_future()
--test_fra()
--test_deposit_rate()
print 'Bootstrapper Pricing Script Loaded Ok'

--fixed_bond_100_6M(39706, 'USD', 'FEDFUND', nil, '15-03-2005:31-08-2010:100.390625', 0.02375)