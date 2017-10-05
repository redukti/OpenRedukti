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

function test_index1()
	local idx = index('USD', 'LIBOR', '1W')
	assert(idx)
	local dt = date(23,10,2016)
	local adjusted = idx:adjust_date(dt, 1)
	assert(adjusted == date(24,10,2016))
	local fixing_dt, value_dt, maturity_dt = idx:date_components(adjusted)
	assert(maturity_dt == date(31,10,2016))
	assert(value_dt == adjusted)
end


test_index1()