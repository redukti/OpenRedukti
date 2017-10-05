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


local utils = assert(require('utils'))
local deposit_rate = utils.deposit_rate

local fixing_service = utils.load_fixings('../testdata/20121211/fixings.csv')

local curve_mapper = redukti.curve_mapper()
local f_eonia_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EONIA'}
local d_eonia_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR' }
local d_euribor_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR'}
local f_euribor_1m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '1M'}
local f_euribor_3m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '3M'}
local f_euribor_6m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '6M'}
local f_euribor_12m_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = '12M'}
local d_euribor_1m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '1M'}
local d_euribor_3m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '3M'}
local d_euribor_6m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '6M'}
local d_euribor_12m_id = redukti.pricing_curve_id { curve_type = 'D', currency = 'EUR', index_family = 'EURIBOR', tenor = '12M'}
local f_euribor_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR'}

-- map request for forward EONIA curve to discount curve
curve_mapper:add_mapping( f_eonia_id, d_eonia_id )
-- map any request for EURIBOR Discount to EONIA discount
curve_mapper:add_mapping( d_euribor_id, d_eonia_id )
-- map euribor 1m to generic
curve_mapper:add_mapping( f_euribor_1m_id, f_euribor_id )
-- map euribor 12m to generic
curve_mapper:add_mapping( f_euribor_1m_id, f_euribor_id )

local id

id = curve_mapper:get_mapping( f_eonia_id )
assert(id == d_eonia_id)

id = curve_mapper:get_mapping( d_euribor_id )
assert(id == d_eonia_id)

id = curve_mapper:get_mapping( d_euribor_3m_id )
assert(id == d_eonia_id)

id = curve_mapper:get_mapping( d_euribor_1m_id )
assert(id == d_eonia_id)

id = curve_mapper:get_mapping( d_euribor_6m_id )
assert(id == d_eonia_id)

local business_date = redukti.date('2012/12/11')
local valuation_context = redukti.valuation_context({ reference_date = business_date, order = 1 }, fixing_service)
assert(debug.getuservalue(valuation_context).fixing_service == fixing_service)

assert(valuation_context)

local flows = deposit_rate(redukti.date(25, 7, 2013), 'EUR', 'EURIBOR', '1Y', 0.00140)
local pricing_cashflows = redukti.prepare_cashflows_for_pricing(valuation_context, curve_mapper, flows)

assert(pricing_cashflows:ok())

local data = redukti.loadcsv { file='../testdata/20121211/zerocurves.csv', conversion='dnnnn', heading=true, default_number = -1.0, fields=true }
local defs = redukti.loadcsv { file='../testdata/20121211/curve_definitions.csv', conversion='sssssssss', heading=true, fields = true }

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
	--print('curve_id = ' .. k)
	--print('currency = ' .. v.currency)
	--print('index_family = ' .. v.index_family)
	--print('tenor = ' .. v.tenor)
	--print('interpolator = ' .. v.interpolator)
	--print('value_type = ' .. v.value_type)
	--print('maturities = ', table.unpack(v.maturities))
	--print('values = ', table.unpack(v.values))

	local zc = redukti.curve(v)
	assert(zc)

	zerocurves[k] = zc
end


local provider = redukti.curve_provider()

for k,v in pairs(zerocurves) do
	local def = curves[k]
	assert(def)

	if def.index_family == 'EONIA' then
		provider:add_mapping(d_eonia_id, v)		
	elseif def.index_family == 'EURIBOR' then
		if def.tenor == nil or def.tenor == '' then
			-- default forward curve
			provider:add_mapping(f_euribor_id, v)
		else 
			local c_id = redukti.pricing_curve_id { curve_type = 'F', currency = 'EUR', index_family = 'EURIBOR', tenor = def.tenor }
			provider:add_mapping(c_id, v)
		end
	end	

end

local zc_eonia = zerocurves['7']
local result = provider:curves()[curve_mapper:get_mapping(f_eonia_id)]
assert(result)
assert(result == zc_eonia)
result = provider:curves()[curve_mapper:get_mapping(d_euribor_6m_id)]
assert(result)
assert(result == zc_eonia)

local pricing_result = redukti.present_value(valuation_context, pricing_cashflows, provider)

assert(pricing_result:ok())

local curve_ids = pricing_result:curve_ids()
assert(curve_ids)

for k,v in pairs(curve_ids) do
	print(k,v)

	local delta = pricing_result:delta(k)
	assert(delta)

	for i = 1,#delta do
		print(i, delta[i])
	end
end

