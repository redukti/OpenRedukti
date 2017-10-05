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

local function compare_sensitivities(a, b)

	local grad1 = a:gradient();
	local grad2 = b:gradient();

	if not grad1 and not grad2 then
		return
	end

	compare_gradient(grad1, grad2)

	local hess1 = a:hessian();
	local hess2 = b:hessian();

	if not hess1 and not hess2 then
		return
	end

	compare_hessian(hess1, hess2)
end

local x : number[] = {
0.01,
0.02,
0.03,
0.04,
0.05,
}

local y : number[] = {
1000000.0,
20004.0,
300000.5,
4000000.0,
900000.0,
}


local interp1 = redukti.interpolator {
	x = x,
	y = y,
	interpolator = 'CubicSplineNatural',
	order = 2
}

local interp2 = redukti.interpolator {
	x = x,
	y = y,
	interpolator = 'LogCubicSplineNatural',
	order = 2
}

-- compare autodiff and numeric sensitivities
compare_sensitivities(
	interp1:interpolate(0.035,1),
	interp1:interpolate(0.035,2))

compare_sensitivities(
	interp2:interpolate(0.035,1),
	interp2:interpolate(0.035,2))

compare_sensitivities(
	interp1:interpolate(0.035,1),
	interp2:interpolate(0.035,1))


local interp3 = redukti.interpolator {
	x = x,
	y = y,
	interpolator = 'Linear',
	order = 2
}

local interp4 = redukti.interpolator {
	x = x,
	y = y,
	interpolator = 'LogLinear',
	order = 2
}

compare_sensitivities(
	interp3:interpolate(0.035,1),
	interp3:interpolate(0.035,2))

compare_sensitivities(
	interp4:interpolate(0.035,1),
	interp4:interpolate(0.035,2))

compare_sensitivities(
	interp3:interpolate(0.035,1),
	interp4:interpolate(0.035,1))

--print(interp3:interpolate(0.035,1))

local interp5 = redukti.interpolator {
	x = x,
	y = y,
	interpolator = 'MonotoneConvex',
	order = 2
}
compare_sensitivities(
	interp5:interpolate(0.035,1),
	interp5:interpolate(0.035,2))
--print(interp5:interpolate(0.035,1))


print 'Ok'