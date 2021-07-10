#include <stdlib.h>
#include <tinyspline.h>
#include "CuTest.h"

#define EPSILON 0.0001

void eval_domain_min(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal *ctrlp = NULL, *result = NULL;
	tsReal min, max;
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			7, 2, 3, TS_CLAMPED, &spline, &status))

		TS_CALL(try, status.code, ts_bspline_control_points(
			&spline, &ctrlp, &status))
		ctrlp[0]  = -1.75f; /* x0 */
		ctrlp[1]  = -1.0f;  /* y0 */
		ctrlp[2]  = -1.5f;  /* x1 */
		ctrlp[3]  = -0.5f;  /* y1 */
		ctrlp[4]  = -1.5f;  /* x2 */
		ctrlp[5]  =  0.0f;  /* y2 */
		ctrlp[6]  = -1.25f; /* x3 */
		ctrlp[7]  =  0.5f;  /* y3 */
		ctrlp[8]  = -0.75f; /* x4 */
		ctrlp[9]  =  0.75f; /* y4 */
		ctrlp[10] =  0.0f;  /* x5 */
		ctrlp[11] =  0.5f;  /* y5 */
		ctrlp[12] =  0.5f;  /* x6 */
		ctrlp[13] =  0.0f;  /* y6 */
		TS_CALL(try, status.code, ts_bspline_set_control_points(
			&spline, ctrlp, &status))

/* ================================= When ================================== */
		ts_bspline_domain(&spline, &min, &max);
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, min, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 1, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 2, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, min, ts_deboornet_knot(&net), EPSILON);
		CuAssertDblEquals(tc, 3, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 4, ts_deboornet_multiplicity(&net),
				  EPSILON);
		CuAssertDblEquals(tc, 0, ts_deboornet_num_insertions(&net),
				  EPSILON);

		TS_CALL(try, status.code, ts_deboornet_result(
			&net, &result, &status))
		CuAssertDblEquals(tc, -1.75f, result[0], EPSILON);
		CuAssertDblEquals(tc, -1.0f, result[1], EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		free(ctrlp);
		ts_deboornet_free(&net);
		free(result);
	TS_END_TRY
}

void eval_domain_max(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal *ctrlp = NULL, *result = NULL;
	tsReal min, max;
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			7, 2, 3, TS_CLAMPED, &spline, &status))

		TS_CALL(try, status.code, ts_bspline_control_points(
			&spline, &ctrlp, &status))
		ctrlp[0]  = -1.75f; /* x0 */
		ctrlp[1]  = -1.0f;  /* y0 */
		ctrlp[2]  = -1.5f;  /* x1 */
		ctrlp[3]  = -0.5f;  /* y1 */
		ctrlp[4]  = -1.5f;  /* x2 */
		ctrlp[5]  =  0.0f;  /* y2 */
		ctrlp[6]  = -1.25f; /* x3 */
		ctrlp[7]  =  0.5f;  /* y3 */
		ctrlp[8]  = -0.75f; /* x4 */
		ctrlp[9]  =  0.75f; /* y4 */
		ctrlp[10] =  0.0f;  /* x5 */
		ctrlp[11] =  0.5f;  /* y5 */
		ctrlp[12] =  0.5f;  /* x6 */
		ctrlp[13] =  0.0f;  /* y6 */
		TS_CALL(try, status.code, ts_bspline_set_control_points(
			&spline, ctrlp, &status))

/* ================================= When ================================== */
		ts_bspline_domain(&spline, &min, &max);
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, max, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 1, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 2, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, max, ts_deboornet_knot(&net), EPSILON);
		CuAssertDblEquals(tc, 10, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 4, ts_deboornet_multiplicity(&net),
				  EPSILON);
		CuAssertDblEquals(tc, 0, ts_deboornet_num_insertions(&net),
				  EPSILON);

		TS_CALL(try, status.code, ts_deboornet_result(
			&net, &result, &status))
		CuAssertDblEquals(tc, 0.5f, result[0], EPSILON);
		CuAssertDblEquals(tc, 0.0f, result[1], EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		free(ctrlp);
		ts_deboornet_free(&net);
		free(result);
	TS_END_TRY
}

void eval_001(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal *ctrlp = NULL, *result = NULL;
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			7, 2, 3, TS_CLAMPED, &spline, &status))

		TS_CALL(try, status.code, ts_bspline_control_points(
			&spline, &ctrlp, &status))
		ctrlp[0]  = -1.75f; /* x0 */
		ctrlp[1]  = -1.0f;  /* y0 */
		ctrlp[2]  = -1.5f;  /* x1 */
		ctrlp[3]  = -0.5f;  /* y1 */
		ctrlp[4]  = -1.5f;  /* x2 */
		ctrlp[5]  =  0.0f;  /* y2 */
		ctrlp[6]  = -1.25f; /* x3 */
		ctrlp[7]  =  0.5f;  /* y3 */
		ctrlp[8]  = -0.75f; /* x4 */
		ctrlp[9]  =  0.75f; /* y4 */
		ctrlp[10] =  0.0f;  /* x5 */
		ctrlp[11] =  0.5f;  /* y5 */
		ctrlp[12] =  0.5f;  /* x6 */
		ctrlp[13] =  0.0f;  /* y6 */
		TS_CALL(try, status.code, ts_bspline_set_control_points(
			&spline, ctrlp, &status))

/* ================================= When ================================== */
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, 0.4f, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 1, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 2, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, 0.4f, ts_deboornet_knot(&net), EPSILON);
		CuAssertDblEquals(tc, 4, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 0, ts_deboornet_multiplicity(&net),
				  EPSILON);
		CuAssertDblEquals(tc, 3, ts_deboornet_num_insertions(&net),
				  EPSILON);

		TS_CALL(try, status.code, ts_deboornet_result(
			&net, &result, &status))
		CuAssertDblEquals(tc, -1.338333f, result[0], EPSILON);
		CuAssertDblEquals(tc,  0.288333f, result[1], EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		free(ctrlp);
		ts_deboornet_free(&net);
		free(result);
	TS_END_TRY
}

void eval_002(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal *ctrlp = NULL, *result = NULL;
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			7, 2, 3, TS_CLAMPED, &spline, &status))

		TS_CALL(try, status.code, ts_bspline_control_points(
			&spline, &ctrlp, &status))
		ctrlp[0]  = -1.75f; /* x0 */
		ctrlp[1]  = -1.0f;  /* y0 */
		ctrlp[2]  = -1.5f;  /* x1 */
		ctrlp[3]  = -0.5f;  /* y1 */
		ctrlp[4]  = -1.5f;  /* x2 */
		ctrlp[5]  =  0.0f;  /* y2 */
		ctrlp[6]  = -1.25f; /* x3 */
		ctrlp[7]  =  0.5f;  /* y3 */
		ctrlp[8]  = -0.75f; /* x4 */
		ctrlp[9]  =  0.75f; /* y4 */
		ctrlp[10] =  0.0f;  /* x5 */
		ctrlp[11] =  0.5f;  /* y5 */
		ctrlp[12] =  0.5f;  /* x6 */
		ctrlp[13] =  0.0f;  /* y6 */
		TS_CALL(try, status.code, ts_bspline_set_control_points(
			&spline, ctrlp, &status))

/* ================================= When ================================== */
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, 0.8f, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 1, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 2, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, 0.8f, ts_deboornet_knot(&net), EPSILON);
		CuAssertDblEquals(tc, 6, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 0, ts_deboornet_multiplicity(&net),
				  EPSILON);
		CuAssertDblEquals(tc, 3, ts_deboornet_num_insertions(&net),
				  EPSILON);

		TS_CALL(try, status.code, ts_deboornet_result(
			&net, &result, &status))
		CuAssertDblEquals(tc, -0.470667f, result[0], EPSILON);
		CuAssertDblEquals(tc,  0.618667f, result[1], EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		free(ctrlp);
		ts_deboornet_free(&net);
		free(result);
	TS_END_TRY
}

void eval_003(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal *ctrlp = NULL, *result = NULL;
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			7, 2, 3, TS_CLAMPED, &spline, &status))

		TS_CALL(try, status.code, ts_bspline_control_points(
			&spline, &ctrlp, &status))
		ctrlp[0]  = -1.75f; /* x0 */
		ctrlp[1]  = -1.0f;  /* y0 */
		ctrlp[2]  = -1.5f;  /* x1 */
		ctrlp[3]  = -0.5f;  /* y1 */
		ctrlp[4]  = -1.5f;  /* x2 */
		ctrlp[5]  =  0.0f;  /* y2 */
		ctrlp[6]  = -1.25f; /* x3 */
		ctrlp[7]  =  0.5f;  /* y3 */
		ctrlp[8]  = -0.75f; /* x4 */
		ctrlp[9]  =  0.75f; /* y4 */
		ctrlp[10] =  0.0f;  /* x5 */
		ctrlp[11] =  0.5f;  /* y5 */
		ctrlp[12] =  0.5f;  /* x6 */
		ctrlp[13] =  0.0f;  /* y6 */
		TS_CALL(try, status.code, ts_bspline_set_control_points(
			&spline, ctrlp, &status))

/* ================================= When ================================== */
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, 0.75f, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 1, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 2, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, 0.75f, ts_deboornet_knot(&net), EPSILON);
		CuAssertDblEquals(tc, 6, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 1, ts_deboornet_multiplicity(&net),
				  EPSILON);
		CuAssertDblEquals(tc, 2, ts_deboornet_num_insertions(&net),
				  EPSILON);

		TS_CALL(try, status.code, ts_deboornet_result(
			&net, &result, &status))
		CuAssertDblEquals(tc, -0.645833f, result[0], EPSILON);
		CuAssertDblEquals(tc,  0.645833f, result[1], EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		free(ctrlp);
		ts_deboornet_free(&net);
		free(result);
	TS_END_TRY
}

void eval_two_points(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal *ctrlp = NULL, *result = NULL;
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			8, 2, 3, TS_BEZIERS, &spline, &status))

		TS_CALL(try, status.code, ts_bspline_control_points(
			&spline, &ctrlp, &status))
		ctrlp[0]  = -1.75f; /* x0 */
		ctrlp[1]  = -1.0f;  /* y0 */
		ctrlp[2]  = -1.5f;  /* x1 */
		ctrlp[3]  = -0.5f;  /* y1 */
		ctrlp[4]  = -1.5f;  /* x2 */
		ctrlp[5]  =  0.0f;  /* y2 */
		ctrlp[6]  = -1.25f; /* x3 */
		ctrlp[7]  =  0.5f;  /* y3 */
		ctrlp[8]  = -0.75f; /* x4 */
		ctrlp[9]  =  0.75f; /* y4 */
		ctrlp[10] =  0.0f;  /* x5 */
		ctrlp[11] =  0.5f;  /* y5 */
		ctrlp[12] =  0.5f;  /* x6 */
		ctrlp[13] =  0.0f;  /* y6 */
		ctrlp[14] = -0.3f;  /* x7 */
		ctrlp[15] = -1.0f;  /* y7 */
		TS_CALL(try, status.code, ts_bspline_set_control_points(
			&spline, ctrlp, &status))

/* ================================= When ================================== */
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, 0.5f, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 2, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 2, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, 0.5f, ts_deboornet_knot(&net), EPSILON);
		CuAssertDblEquals(tc, 7, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 4, ts_deboornet_multiplicity(&net),
				  EPSILON);
		CuAssertDblEquals(tc, 0, ts_deboornet_num_insertions(&net),
				  EPSILON);

		TS_CALL(try, status.code, ts_deboornet_result(
			&net, &result, &status))
		CuAssertDblEquals(tc, -1.25, result[0], EPSILON);
		CuAssertDblEquals(tc,  0.5f, result[1], EPSILON);
		CuAssertDblEquals(tc, -0.75f, result[2], EPSILON);
		CuAssertDblEquals(tc,  0.75f, result[3], EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		free(ctrlp);
		ts_deboornet_free(&net);
		free(result);
	TS_END_TRY
}

void eval_undefined_knot(CuTest *tc)
{
	const tsReal step = (TS_DOMAIN_DEFAULT_MAX - TS_DOMAIN_DEFAULT_MIN) /
		(tsReal) 10.f; /* num(control_points) + degree */
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsStatus status;

/* ================================= Given ================================= */
	CuAssertIntEquals(tc, TS_SUCCESS, ts_bspline_new(
		7, 3, 3, TS_OPENED, &spline, &status));

/* =============================== When/Then =============================== */
	CuAssertIntEquals(tc, TS_U_UNDEFINED, ts_bspline_eval(
		&spline, TS_DOMAIN_DEFAULT_MIN, &net, &status));
	CuAssertIntEquals(tc, TS_U_UNDEFINED, ts_bspline_eval(
		&spline, step, &net, &status));
	CuAssertIntEquals(tc, TS_U_UNDEFINED, ts_bspline_eval(
		&spline, step * 2, &net, &status));
	CuAssertIntEquals(tc, TS_U_UNDEFINED, ts_bspline_eval(
		&spline, step * 8, &net, &status));
	CuAssertIntEquals(tc, TS_U_UNDEFINED, ts_bspline_eval(
		&spline, step * 9, &net, &status));
	CuAssertIntEquals(tc, TS_U_UNDEFINED, ts_bspline_eval(
		&spline, TS_DOMAIN_DEFAULT_MAX, &net, &status));
}

void eval_near_miss_knot(CuTest *tc)
{
	tsBSpline spline = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsStatus status;

	TS_TRY(try, status.code, &status)
/* ================================= Given ================================= */
		TS_CALL(try, status.code, ts_bspline_new(
			3, 4, 2, TS_CLAMPED, &spline, &status))
		CuAssertIntEquals(tc, 6, ts_bspline_num_knots(&spline));

		TS_CALL(try, status.code, ts_bspline_set_knots_varargs(
			&spline, &status,
			0.1, 0.2, 0.3, 0.3, 0.4, 0.6))

/* ================================= When ================================== */
		TS_CALL(try, status.code, ts_bspline_eval(
			&spline, 0.2999999f, &net, &status))

/* ================================= Then ================================== */
		CuAssertIntEquals(tc, 1, ts_deboornet_num_result(&net));
		CuAssertIntEquals(tc, 4, ts_deboornet_dimension(&net));
		CuAssertDblEquals(tc, 0.3, ts_deboornet_knot(&net), 0);
		CuAssertDblEquals(tc, 3, ts_deboornet_index(&net), EPSILON);
		CuAssertDblEquals(tc, 2, ts_deboornet_multiplicity(&net),
			EPSILON);
		CuAssertDblEquals(tc, 0, ts_deboornet_num_insertions(&net),
			EPSILON);
	TS_CATCH(status.code)
		CuFail(tc, status.message);
	TS_FINALLY
		ts_bspline_free(&spline);
		ts_deboornet_free(&net);
	TS_END_TRY
}

CuSuite* get_eval_suite()
{
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, eval_domain_min);
	SUITE_ADD_TEST(suite, eval_domain_max);
	SUITE_ADD_TEST(suite, eval_001);
	SUITE_ADD_TEST(suite, eval_002);
	SUITE_ADD_TEST(suite, eval_003);
	SUITE_ADD_TEST(suite, eval_two_points);
	SUITE_ADD_TEST(suite, eval_undefined_knot);
	SUITE_ADD_TEST(suite, eval_near_miss_knot);
	return suite;
}