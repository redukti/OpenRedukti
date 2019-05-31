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


-- This particular example is taken from Examples/Bond.cpp
-- in Quantlib ; our curve doesn't quite match what Quantlib generates

local business_date = redukti.date('2008/09/18')
local utils = assert(require('utils'))
local result, errmsg = utils.build_curves(business_date, 
	'../testdata/bond/bond_curve_definitions.txt', 
	'../testdata/bond/bond_curve_inputs.txt',
	{1})
if not result:ok() then print(v) end
assert(result:ok())

print()
local curves = result:curves()
for k,v in pairs(curves) do
	local maturities, rates, factors = v:values()
	for i = 1,#maturities do 
		print('Maturity ' .. maturities[i] .. ' rate ' .. rates[i] .. ' discount factor ' .. factors[i])
	end
end

local expected_maturities = {}
local expected_rates = {}
local expected_factors = {}

expected_maturities[1] = {
	39800,
	39890,
	40074,
	40421,
	40786,
	41517,
	43327,
	50540
}
expected_rates[1] = {
	0.009588529859942798,
	0.014448118073675386,
	0.019214218923804424,
	0.022333814814160188,
	0.025213394023920782,
	0.03041127799766808,
	0.039079008321177056,
	0.04645988031116752
}
expected_factors[1] = {
	0.9976122901461,
	0.99286092194617,
	0.9809691975672,
	0.95736913966835,
	0.92830328689121,
	0.86015732166091,
	0.67884396543265,
	0.25191861109561
}

-- Above doesn't match what Quantlib generates:
--Curve node #0 date: 39709 value: 1
--Curve node #1 date: 39800 value: 0.997612
--Curve node #2 date: 39890 value: 0.992861
--Curve node #3 date: 40074 value: 0.980969
--Curve node #4 date: 40421 value: 0.958714
--Curve node #5 date: 40786 value: 0.930809
--Curve node #6 date: 41517 value: 0.861715
--Curve node #7 date: 43327 value: 0.681765
--Curve node #8 date: 50540 value: 0.260189

local function compare(a,b)
	for i = 1,#a do
		assert(math.abs(a[i]-b[i]) < 1e-10)
	end
end

local curves = result:curves()
for k,v in pairs(curves) do
	print('Checking ' .. k .. ' ' .. tostring(v))

	local maturities, rates, factors = v:values()

	compare(maturities, expected_maturities[k])
	compare(rates, expected_rates[k])
	compare(factors, expected_factors[k])
end
