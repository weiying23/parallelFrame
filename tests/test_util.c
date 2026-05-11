#include "mythread/mythread_util.h"
#include "mythread/mythread_types.h"
#include <stdio.h>
#include <string.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

/* TC-U01: ntdelay does not crash */
static void test_ntdelay(void) {
  TEST("ntdelay-no-crash");
  ntdelay(1);
  ntdelay(100);
  ntdelay(0);
  ntdelay(-1);
  ntdelay(10000);
  PASS();
}

/* TC-U02: ntdelay_ does not crash (data dependency loop) */
static void test_ntdelay_(void) {
  TEST("ntdelay_-no-crash");
  ntdelay_(1);
  ntdelay_(100);
  ntdelay_(0);
  PASS();
}

/* TC-U03: opentf does not crash */
static void test_opentf(void) {
  TEST("opentf-no-crash");
  /* without THLOG defined, opentf is a no-op */
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  t.fo = NULL;
  ti = &t;
  opentf();
  CHK(t.fo == NULL, "opentf without THLOG should not set fo");
  ti = saved_ti;
  PASS();
}

/* TC-U04: HNEW macro works */
static void test_hnew_macro(void) {
  TEST("hnew-macro");
  THREADINFO *t = HNEW(THREADINFO);
  CHK(t != NULL, "HNEW returned NULL");
  memset(t, 0, sizeof(*t));
  free(t);
  PASS();
}

int main(void) {
  printf("\n=== mythread util tests ===\n\n");

  test_ntdelay();
  test_ntdelay_();
  test_opentf();
  test_hnew_macro();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
