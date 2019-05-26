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

local function table_print (tt, indent, done)
  done = done or {}
  indent = indent or 0
  if type(tt) == "table" then
    for key, value in pairs (tt) do
      io.write(string.rep (" ", indent)) -- indent it
      if type (value) == "table" and not done [value] then
        done [value] = true
        io.write(string.format("[%s] => table\n", tostring (key)));
        io.write(string.rep (" ", indent+4)) -- indent it
        io.write("(\n");
        table_print (value, indent + 7, done)
        io.write(string.rep (" ", indent+4)) -- indent it
        io.write(")\n");
      else
        io.write(string.format("[%s] => %s\n",
          tostring (key), tostring(value)))
      end
    end
  else
    io.write(tt .. "\n")
  end
end

local function load_fixings(filename)
	local fixings = redukti.loadcsv { file=filename, conversion='ssdn', heading=true, fields=true }

	assert(fixings)

	local indices = {}
	local fixings_by_index = {}

	for i=1,#fixings do
		local f = fixings[i]
		local indexkey = f.index .. f.tenor
		local index = indices[indexkey]
		if not index then
			index = redukti.index(f.index, f.tenor)
			assert(index)
			indices[indexkey] = index
		end
		local data = fixings_by_index[index]
		if not data then
			data = {}
			fixings_by_index[index] = data
		end
		data[f.date] = f.fixing
	end

	return redukti.fixing_service(fixings_by_index), fixings_by_index
end

local function deposit_rate(business_date: integer, ccy, index_family, tenor, rate: number)
	local day_count = 'ACT/365.FIXED' 
	local start_delay = 0
	if (ccy ~= 'GBP') then
		if (tenor == '1D') then
			start_delay = 0
		elseif (tenor == '2D') then
			start_delay = 1
		else
			start_delay = 2
		end
		day_count = 'ACT/360'
	end
	local index = redukti.index(ccy, index_family, tenor)
	if not index then
		return nil
	end
	local notional: number = 1000000.0
	local start_date: integer = index:adjust_date(business_date, start_delay)
	local maturity_date: integer = index:maturity(start_date)
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
	return flows
end

local function equals(a: number, b: number)
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

-- search for value in list
local function find(list, value)
	for i = 1,#list do
		if list[i] == value then
			return true
		end
	end
	return false
end

-- only include definitions that have matching ids in keeplist
local function filter_definitions(defs, keeplist)
	local newdefs = {}
	for i = 1,#defs do
		if find(keeplist, defs[i].curve_id) then 
			table.insert(newdefs, defs[i])
		end
	end
	print(table.unpack(newdefs))
	return newdefs
end

local function build_curves(business_date, defs_file, par_rates_file, curveid_list)
	local defs = redukti.loadcsv { file=defs_file, conversion='sisssssss', heading=true, fields=true }
	local par_rates = redukti.loadcsv { file=par_rates_file, conversion='isssnii', heading=true, fields=true }
	local result, msg = redukti.build_curves( business_date, filter_definitions(defs, curveid_list), par_rates) -- , '/OpenRedukti/scripts/pricing.lua' )
	if not result:ok() then
		print(msg)
	end
	return result
end

return {
	load_fixings=load_fixings,
	deposit_rate=deposit_rate,
	equals = equals,
	compare_gradient = compare_gradient,
	compare_hessian = compare_hessian,
	build_curves = build_curves,
	table_print = table_print
}