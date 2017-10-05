/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#ifndef _REDUKTI_STATUS_H
#define _REDUKTI_STATUS_H

#include <enums.pb.h>
#include <stddef.h>

namespace redukti
{

typedef ResponseSubCode StatusCode;

extern const char *error_message(StatusCode code);
extern const char *error_message(char *buf, size_t buflen, StatusCode status_code, const char *format, ...);
extern int test_status();
}

#endif
