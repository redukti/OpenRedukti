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

function test_schedule1()

	local t = {}
	t.effective_date = date{1,1,2016}
	t.termination_date = date{1,1,2017}
	t.payment_frequency = "3M"
	t.payment_business_centers = "GBLO,USNY"
	t.payment_day_convention = "MODFOLLOWING"

	local starts, ends, pays = redukti.schedule(t)

	assert(starts ~= nil and #starts == 4)
	assert(ends ~= nil and #ends == 4)
	assert(pays ~= nil and #pays == 4)

	for k,v in pairs(starts) do
		print(k, date_parts(v))
	end
end


test_schedule1()