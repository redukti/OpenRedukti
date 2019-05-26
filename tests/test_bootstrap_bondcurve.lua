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

local business_date = redukti.date('2008/09/15')
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
	39799,
	39889,
	40073,
	40421,
	40786,
	41517,
	43327,
	50540
}
expected_rates[1] = {
	0.0095860765553337,
	0.014592250678254,
	0.019424542573528,
	0.022192046052709,
	0.025084854386565,
	0.030323972407618,
	0.039050572523723,
	0.04645600301827
}
expected_factors[1] = {
	0.99756050064685,
	0.99271058289641,
	0.98065851546356,
	0.95745927620267,
	0.92846399296358,
	0.86031493716342,
	0.67881742407796,
	0.2518514143512
}

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
