local usny_calendar = redukti.calendar('USNY')
local gblo_calendar = redukti.calendar('GBLO')

local usnygblo_calendar = redukti.calendar('USNY,GBLO')
local gblousny_calendar = redukti.calendar('GBLO,USNY')

assert(tostring(usnygblo_calendar) == tostring(gblousny_calendar)) -- check that both map to the same calendar

local gblousny2_calendar = redukti.calendar('GBLO,USNY', 'JOIN_BUSINESS_DAYS')
print(gblousny2_calendar)
