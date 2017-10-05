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
local fs, fixings_by_index = utils.load_fixings('../testdata/20121211/fixings.csv')
assert(fs)

-- search for all the fixing we loaded to make sure
-- thay they are there
for k,v in pairs(fixings_by_index) do
	for d,f in pairs(v) do
		local temp = fs:fixing(k,d)
		assert(temp and temp == f)
	end

	-- now search something that should not be there
	assert(not fs:fixing(k, redukti.date('1970/01/01')))
	assert(not fs:fixing(k, redukti.date('2080/01/01')))
	assert(not fs:fixing(k, redukti.date('2010/12/25')))
end


