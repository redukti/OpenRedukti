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
local equals, compare_hessian, compare_gradient = utils.equals, utils.compare_hessian, utils.compare_gradient
local adouble1 = redukti.adouble1
local adouble2 = redukti.adouble2
local x, y, z, t
local g, h

x = adouble1(2.0)
assert(x() == 2.0)
g = x:gradient()
assert(g and #g == 1 and g[1] == 1.0)
h = x:hessian()
assert(not h)	-- no hessian

x,y = adouble2 {4.2,3.1}
assert(x() == 4.2)
assert(y() == 3.1)
g = x:gradient()
assert(g and #g == 2 and g[1] == 1.0 and g[2] == 0.0)
h = x:hessian()
assert(h, {{0.0, 0.0}, {0.0, 0.0}})
g = y:gradient()
assert(g and #g == 2 and g[1] == 0.0 and g[2] == 1.0)
h = y:hessian()
assert(h, {{0.0, 0.0}, {0.0, 0.0}})

z = x-y
assert(z() == (4.2-3.1))
g = z:gradient()
assert(g and #g == 2 and g[1] == 1.0 and g[2] == -1.0)
h = z:hessian()
assert(h, {{0.0, 0.0}, {0.0, 0.0}})

x = adouble1(2.0)
y = 3*x+5;
assert(y() == 11.0) 
assert(y(1) == 3.0) -- first derivative

x, y, z = adouble2 {5.0, 3.0, 6.0}

assert(x() == 5.0)
assert(y() == 3.0)
assert(z() == 6.0)

t = x+y+z
assert(t() == 14.0)
assert(t(1) == 1.0)
assert(t(2) == 1.0)
assert(t(3) == 1.0)

t = x*y*z
assert(t() == 90.0)
assert(t(1) == 3.0*6.0) -- d/dx = yz
assert(t(2) == 5.0*6.0) -- d/dy = xz
assert(t(3) == 5.0*3.0) -- d/dz = xy
h = t:hessian()
compare_hessian(h, {
	{0.0, 6.0, 3.0},
	{6.0, 0.0, 5.0},
	{3.0, 5.0, 0.0}
	}) 

x,y = adouble2(5.0, 3.0)
z = (x*x)/(y*y*y)

assert(equals(z(), 0.92592592592593))
compare_gradient(z:gradient(), 
	{ 0.37037037037037, -0.92592592592593 })
compare_hessian(z:hessian(), 
	{{ 0.074074074074074, -0.37037037037037 },
    { -0.37037037037037,  1.2345679012346 }})

z = x:pow(2) / y:pow(3)
assert(equals(z(), 0.92592592592593))
compare_gradient(z:gradient(), 
	{ 0.37037037037037, -0.92592592592593 })
compare_hessian(z:hessian(), 
	{{ 0.074074074074074, -0.37037037037037 },
    { -0.37037037037037,  1.2345679012346 }})

