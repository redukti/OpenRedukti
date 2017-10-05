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


local data = redukti.loadcsv { file='../testdata/20121211/zerocurves.csv', conversion='dnnnn', heading=true, default_number = -1.0, fields=true }

local defs = redukti.loadcsv { file='../testdata/20121211/curve_definitions.csv', conversion='sssssssss', heading=true, fields = true }

local function equals(a,b)
	return math.abs(a-b) < 1e-12
end

local function compare_hessian(a, b)
	assert(a)
	assert(b)
	assert(#a == #b)
	for i = 1,#a do
		local row1 = a[i]
		local row2 = b[i]
		assert(#row1 == #row2)
		for j = 1,#row1 do
			assert(equals(row1[j],  row2[j]))
		end
	end
end

local function compare_gradient(a, b)
	assert(a)
	assert(b)
	assert(#a == #b)
	for i = 1,#a do
		assert(equals(a[i], b[i]))
	end
end	

local function defsbyid(defs)
	-- convert from array of definitions
	-- to table of definitions keyed by curve id
	local newdefs = {}
	for i = 1,#defs do
		newdefs[defs[i].curve_id] = defs[i]
	end
	return newdefs
end

defs = defsbyid(defs)

local curves = {}
local daycount = redukti.dayfraction('ACT/365.FIXED')
assert(daycount)
local business_date = redukti.date('2012/12/11')
assert(business_date)

for i = 1,#data do
	local maturity = data[i].maturity

	for k,v in pairs(data[i]) do
		if k ~= 'maturity' then
			if v ~= -1.0 then
				local defn = defs[k]
				assert(defn)

				local c = curves[k]
				if not c then
					c = { maturities = {}, 
						values = {}, 
						currency = defn.currency, 
						index_family = defn.index_family,
						interpolator = defn.interpolator, 
						value_type = defn.interpolated_on, -- we may need to convert, see below
						tenor = defn.curve_tenor,
						reference_date = business_date,
						order = 1
					}
					curves[k] = c
				end
				assert(maturity)

				table.insert(c.maturities, maturity)
				if defn.interpolated_on == 'DiscountFactor' then
					local t = daycount:fraction(business_date, maturity)
					local df = math.exp(-v*t)	-- convert to discountfactor
					table.insert(c.values, df)
				elseif defn.interpolated_on == 'ZeroRate' then
					table.insert(c.values, v)
				else
					error('unsupported value type ' .. defn.interpolated_on)
				end
			end
		end
	end
end

local zerocurves = {}

for k,v in pairs(curves) do
	print('curve_id = ' .. k)
	print('currency = ' .. v.currency)
	print('index_family = ' .. v.index_family)
	print('tenor = ' .. v.tenor)
	print('interpolator = ' .. v.interpolator)
	print('value_type = ' .. v.value_type)
	print('maturities = ', table.unpack(v.maturities))
	print('values = ', table.unpack(v.values))

	local zc = redukti.curve(v)
	assert(zc)

	local d1 = business_date+1;
	local d2 = business_date+10750

	--print(d1,d2)

	for i = 1,100 do

		local d = math.random(d1, d2)

		-- print('Testing date ', d)

		local ad1 = zc:sensitivities(d, 1)
		local ad2 = zc:sensitivities(d, 2)

		compare_gradient(ad1:gradient(), ad2:gradient())
	end

end
