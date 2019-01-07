-- Exercise the api for dates
-- Also acts as a demo

-- Various ways to make a date
local date1 = redukti.date(28,2,2017)
local date2 = redukti.date '2017/02/28'
local date3 = redukti.date '28/02/2017'
local date4 = redukti.date {28,2,2017}

assert(date1 == date2)
assert(date2 == date3)
assert(date3 == date4)

-- Extract date parts
local d,m,y = redukti.date_parts(date1)
assert(d == 28)
assert(m == 2)
assert(y == 2017)

-- The date making functions work nicely with out of range values
date1 = redukti.date(29,2,2017) -- this should change to 1 March 2017 as there is no 29th Feb
d,m,y = redukti.date_parts(date1)

assert(d == 1)
assert(m == 3)
assert(y == 2017)

date1 = redukti.date(0,3,2017)
d,m,y = redukti.date_parts(date1)

assert(d == 28)
assert(m == 2)
assert(y == 2017)

date1 = redukti.addperiod(date1, '-1D') -- move back one day
assert(date1 == date2-1)

assert(not pcall(redukti.date, 'Invalid'))

print 'Date Tests Ok'