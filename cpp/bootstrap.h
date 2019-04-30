/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017-2019 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#ifndef _REDUKTI_BOOTSTRAP_H
#define _REDUKTI_BOOTSTRAP_H

#include <bootstrap.pb.h>
#include <curve.pb.h>
#include <enums.pb.h>

#include <memory>
#include <string>

namespace redukti
{

class CurveBuilderService
{
	public:
	virtual ~CurveBuilderService() {}
	virtual BootstrapCurvesReply *handle_bootstrap_request(google::protobuf::Arena *arena,
							       const BootstrapCurvesRequest *request) = 0;
};

std::unique_ptr<CurveBuilderService> get_curve_builder_service(std::string script);

extern int test_bootstrap();

} // namespace redukti

#endif