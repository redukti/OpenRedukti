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

#include <status.h>

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

namespace redukti
{

const char *error_message(StatusCode code)
{
	switch (code) {
	case StatusCode::kOk:
		return "";
	case StatusCode::kBadArgument:
		return "ERR0001: Bad argument";
	case StatusCode::kError:
		return "ERR0002: Error servicing the request";
	case StatusCode::kNotImplemented:
		return "ERR0003: Feature not yet implemented";
	case StatusCode::kNotAvailable:
		return "ERR0004: Requested service is not available";
	case StatusCode::kInternalError:
		return "ERR0005: Internal error (contact support)";
	case StatusCode::kShuttingDown:
		return "ERR0006: Service is shutting down";
	case StatusCode::kBTS_LuaInitFailure:
		return "BTS1201: Failed to initialize Lua";
	case StatusCode::kBTS_BadInputCurves:
		return "BTS1202: Error processing input par rates";
	case StatusCode::kBTS_InstrumentError:
		return "BTS1203: Error setting up par instruments";
	case StatusCode::kBTS_SolverFailure:
		return "BTS1204: Solver failure";
	case StatusCode::kBTS_BadBusinessDate:
		return "BTS1205: Bad business date";
	case StatusCode::kBTS_BadInterpolatedOn:
		return "BTS1206: Bad interpolated on value";
	case StatusCode::kBTS_MissingFixedTenors:
		return "BTS1207: Curve Definition is missing tenor list for tenor based "
		       "maturity generation";
	case StatusCode::kBTS_BadMaturitiesSortOrder:
		return "BTS1208: The curve maturities are not in correct order or have "
		       "dates < business date";

	case StatusCode::kBTS_CurveDefinitionNotFound: //  = 1209;
		return "BTS1209: Definition for curve not found";
	case StatusCode::kBTS_CashflowDefinitionFailed: // = 1210;
		return "BTS1210: Failed to define cashflows for the instrument";
	case StatusCode::kBTS_LuaScriptFailure: // = 1211;
		return "BTS1211: Lua script failed";
	case StatusCode::kBTS_LuaPricingRoutineNotFound: // = 1212;
		return "BTS1212: Lua pricing routine not found";
	case StatusCode::kBTS_DuplicateCurve: // = 1213;
		return "BTS1213: Curve has same definition id as another curve";
	case StatusCode::kBTS_ForwardCurveReferenceNotFound: // = 1214;
		return "BTS1214: Reference to forward curve not found";
	case StatusCode::kBTS_DiscountCurveReferenceNotFound: // = 1215;
		return "BTS1215: Reference to discount curve not found";
	case StatusCode::kBTS_CashflowGenerationFailed: // = 1216;
		return "BTS1216: Cashflow generation failed";
	case StatusCode::kBTS_MissingCurveDefinitions:
		return "BTS1217: Request is missing curve definitions";
	case StatusCode::kBTS_MissingParCurves:
		return "BTS1218: Request is missing PAR curves";
	case StatusCode::kBTS_DefinitonsAndCurvesMismatch:
		return "BTS1219: Curve Definitions and PAR Curves do not match";
	case StatusCode::kBTS_BadMaturityGenerationRule:
		return "BTS1220: Maturity generation rule is not compatible with curve type";

	case StatusCode::kVAL_BadBusinessDate:
		return "VAL1231: Bad business date";
	case StatusCode::kVAL_BadMaturityDate:
		return "VAL1232: Bad maturity date";
	case StatusCode::kVAL_BadZeroRate:
		return "VAL1233: Zero rate < -0.05 or > 0.2 is probably erroneous";
	case StatusCode::kVAL_CurveDefinitionNotRegistered:
		return "VAL1234: Curve Definition for Zero Curve has not been found; was "
		       "it registered?";
	case StatusCode::kVAL_MismatchedMaturitiesAndValues:
		return "VAL1235: Mismatched number of matuties and values in provided "
		       "curve";
	case StatusCode::kVAL_InsufficientMaturities:
		return "VAL1236: Insufficient maturities in supplied curve; at least 4 "
		       "expected";
	case StatusCode::kVAL_BadFixingDate:
		return "VAL1237: Bad fixing date";
	case StatusCode::kVAL_BadFixingValue:
		return "VAL1238: Fixing value < -0.1 or > 0.2 is probably erroneous";
	case StatusCode::kVAL_UnsupportedIndexOrTenor:
		return "VAL1239: Unsupported Index or Index Tenor";
	case StatusCode::kVAL_MissingCashflow:
		return "VAL1240: Missing cashflow";
	case StatusCode::kVAL_MissingPricingContext:
		return "VAL1241: Missing Pricing Context";
	case StatusCode::kVAL_FixingDataUnavailable:
		return "VAL1242: Fixings are not available";
	case StatusCode::kVAL_FixingNotFound:
		return "VAL1243: Fixing was not found";
	case StatusCode::kVAL_MappedCurveNotFound:
		return "VAL1244: Mapped curve not found";
	case StatusCode::kVAL_UnsupportedCashflowType:
		return "VAL1245: Unsupported cashflow type";
	case StatusCode::kVAL_CurveGroupNotFound:
		return "VAL1246: Specified Curve Group not available, have the curves "
		       "for this group been published?";
	case StatusCode::kVAL_CashflowConversionFailed:
		return "VAL1247: Cashflow conversion failed";
	case StatusCode::kVAL_CurveMappingsForCurveGroupNotFound:
		return "VAL1248: Curve mappings not found for specified curve group";
	case StatusCode::kVAL_DetectedTooManyRiskFactorsInCashflow:
		return "VAL1249: Detected more than 3 risk factors in a cashflow; curve "
		       "mapping may be incorrectly specified";
	case StatusCode::kVAL_MissingParSensitivitiesForCurve:
		return "VAL1250: PAR sensitivities are missing for the zero curve";
	case StatusCode::kVAL_ParSensitivitiesForCurveMismatchedMaturities:
		return "VAL1251: PAR sensitivities do not have same maturities range as "
		       "zero curve";
	case StatusCode::kVAL_ParSensitivitiesForCurveInstrumentsOutOfRange:
		return "VAL1252: PAR sensitivities has too many or too few instruments";

	case StatusCode::kCFG_ScheduleGenerationFailed:
		return "CFG1260: Schedule generation failed";
	case StatusCode::kCFG_MissingTrade:
		return "CFG1261: Invalid input - missing trade";
	case StatusCode::kCFG_UnsupportedLegType:
		return "CFG1262: Unsupported leg type";
	case StatusCode::kCFG_BadInput:
		return "CFG1263: Trade failed validation";
	case StatusCode::kCFG_UnsupportedTradeType:
		return "CFG1264: Unsupported trade type";
	case StatusCode::kCFG_SteppedCashflowsNotImplementedYet:
		return "CFG1265: Stepped cashflows not implemented yet";
	case StatusCode::kCFG_IndexConfigurationMissing:
		return "CFG1266: Index configuration not found";

	case StatusCode::kSCH_FirstRegularStartDateInconsistent:
		return "SCH1300: First regular start date is not consistent with the Roll Convention";
	case StatusCode::kSCH_LastRegularEndDateInconsistent:
		return "SCH1301: Last regular period end date is not consistent with the Roll Convention";
	case StatusCode::kSCH_EffectiveDateRequired:
		return "SCH1302: Effective date is required";
	case StatusCode::kSCH_TerminateDateOrTermRequired:
		return "SCH1303: One of termination date or term is required";
	case StatusCode::kSCH_TerminationDateRequired:
		return "SCH1304: Termination date is required";
	case StatusCode::kSCH_CalculationFrequencyRequired:
		return "SCH1305: Calculation frequency is required";
	case StatusCode::kSCH_IncompatiblePaymentAndCalculationFrequencies:
		return "SCH1306: Incompatible payment and calculation frequencies";
	case StatusCode::kSCH_PaymentFrequencyLessThanCalculationFrequency:
		return "SCH1307: Payment frequency must be greater than calculaton frequency";
	case StatusCode::kSCH_PaymentFrequencyNotMultipleOfCalculationFrequency:
		return "SCH1308: Payment frequency must be a multiple of the calculation frequency";

	case StatusCode::kCAL_RegisterCalendarFailed:
		return "CAL1401: Failed to register calendar: note calendars can only be registered before calendar "
		       "service is used";

	case StatusCode::kIDX_ValidFixingCalendarRequired:
		return "IDX1410: Failed to register index definition because definition lacks valid fixing calendar";
	case StatusCode::kIDX_ValidDayCountFractionRequired:
		return "IDX1411: Failed to register index definition because definition lacks valid day count fraction";

	default:
		return "ERR0000: unexpected error";
	}
}

const char *error_message(char *buf, size_t buflen, StatusCode status_code, const char *format, ...)
{
	const char *status_msg = error_message(status_code);
	int n = snprintf(buf, buflen, "%s", status_msg);
	if (n < 0 || n >= buflen) {
		return buf;
	}
	char *buf2 = buf + n;
	buflen -= n;
	va_list args;
	va_start(args, format);
	vsnprintf(buf2, buflen, format, args);
	va_end(args);
	return buf;
}

int test_status()
{
	int failure_count = 0;
	char buf[128];
	const char *msg = error_message(buf, sizeof buf, StatusCode::kCFG_IndexConfigurationMissing, ": %d, %d", 5, 6);
	if (strcmp(msg, "CFG1266: Index configuration not found: 5, 6") != 0)
		failure_count++;
	if (failure_count == 0)
		printf("Test Status OK\n");
	else
		printf("Test Status FAILED\n");
	return failure_count;
}
} // namespace redukti
