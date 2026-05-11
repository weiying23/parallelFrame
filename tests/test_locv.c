#include "mythread/mythread.h"
#include <stdio.h>
#include <string.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

/* TC-L01: SetLocV/GetLocV for typ=2 (md) with low indices */
static void test_locv_md_low(void) {
  TEST("locv-md-low");
  short saved_nlocv = md.nlocv;
  void *saved_locv0 = md.locv[0];
  void *saved_locv1 = md.locv[1];
  void *saved_locv7 = md.locv[7];

  md.nlocv = 8;  /* ensure enough space */

  int dummy0 = 42, dummy1 = 99;
  SetLocV(2, 0, &dummy0);
  SetLocV(2, 1, &dummy1);

  CHK(GetLocV(2, 0, NULL) == &dummy0, "md locv[0] mismatch");
  CHK(GetLocV(2, 1, NULL) == &dummy1, "md locv[1] mismatch");
  CHK(*(int*)GetLocV(2, 0, NULL) == 42, "md locv[0] value mismatch");

  /* restore */
  md.locv[0] = saved_locv0;
  md.locv[1] = saved_locv1;
  md.locv[7] = saved_locv7;
  md.nlocv = saved_nlocv;
  PASS();
}

/* TC-L02: GetLocV returns NULL for unset indices */
static void test_locv_get_unset(void) {
  TEST("locv-get-unset");
  short saved_nlocv = md.nlocv;
  void *saved_locv3 = md.locv[3];

  md.nlocv = 3;  /* only indices 0-2 valid */
  md.locv[3] = NULL;

  CHK(GetLocV(2, 3, NULL) == NULL, "unset md locv[3] should be NULL");

  md.locv[3] = saved_locv3;
  md.nlocv = saved_nlocv;
  PASS();
}

/* TC-L03: SetLocV/GetLocV for typ=1 (gi) with valid gi */
static void test_locv_gi(void) {
  TEST("locv-gi");
  threadGroup *saved_gi = gi;
  threadGroup g;
  memset(&g, 0, sizeof(g));
  g.nlocv = 8;
  gi = &g;

  int value = 123;
  SetLocV(1, 0, &value);
  CHK(GetLocV(1, 0, NULL) == &value, "gi locv[0] mismatch");

  gi = saved_gi;
  PASS();
}

/* TC-L04: SetLocV/GetLocV for typ=0 (ti) with valid ti */
static void test_locv_ti(void) {
  TEST("locv-ti");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  t.nlocv = 8;
  ti = &t;

  int value = 456;
  SetLocV(0, 0, &value);
  CHK(GetLocV(0, 0, NULL) == &value, "ti locv[0] mismatch");

  ti = saved_ti;
  PASS();
}

/* TC-L05: SetLocV/GetLocV return NULL for invalid typ */
static void test_locv_invalid_typ(void) {
  TEST("locv-invalid-typ");
  int dummy = 1;
  /* typ=3 is invalid */
  SetLocV(3, 0, &dummy);  /* should not crash */
  CHK(GetLocV(3, 0, NULL) == NULL, "invalid typ should return NULL");
  PASS();
}

/* TC-L06: SetLocV/GetLocV reject negative index */
static void test_locv_negative_index(void) {
  TEST("locv-negative-index");
  int dummy = 1;
  SetLocV(2, -1, &dummy);  /* should not crash */
  CHK(GetLocV(2, -1, NULL) == NULL, "negative index should return NULL");
  PASS();
}

/* TC-L07: SetLocV/GetLocV reject index > 100 */
static void test_locv_high_index(void) {
  TEST("locv-high-index");
  int dummy = 1;
  SetLocV(2, 101, &dummy);  /* should not crash */
  CHK(GetLocV(2, 101, NULL) == NULL, "index > 100 should return NULL");
  PASS();
}

/* TC-L08: SetLocV with NULL ti/gi returns safely */
static void test_locv_null_context(void) {
  TEST("locv-null-context");
  THREADINFO *saved_ti = ti;
  threadGroup *saved_gi = gi;
  ti = NULL;
  gi = NULL;

  int value = 5;
  SetLocV(0, 0, &value);  /* ti is NULL, should return safely */
  CHK(GetLocV(0, 0, NULL) == NULL, "NULL ti should return NULL");

  SetLocV(1, 0, &value);  /* gi is NULL, should return safely */
  CHK(GetLocV(1, 0, NULL) == NULL, "NULL gi should return NULL");

  ti = saved_ti;
  gi = saved_gi;
  PASS();
}

/* TC-L09: set then get same value across all typ levels */
static void test_locv_set_get_roundtrip(void) {
  TEST("locv-set-get-roundtrip");
  short saved_md_nlocv = md.nlocv;
  void *saved_locv[8];
  memcpy(saved_locv, md.locv, sizeof(saved_locv));

  md.nlocv = 8;
  int v = 777;
  int v_check = 0;
  CHK(GetLocV(2, 2, NULL) != &v, "pre-condition: should not already point to v");
  SetLocV(2, 2, &v);
  CHK(GetLocV(2, 2, NULL) == &v, "roundtrip failed");
  CHK(*(int*)GetLocV(2, 2, NULL) == 777, "roundtrip value mismatch");

  memcpy(md.locv, saved_locv, sizeof(saved_locv));
  md.nlocv = saved_md_nlocv;
  PASS();
}

int main(void) {
  printf("\n=== mythread locv tests ===\n\n");

  test_locv_md_low();
  test_locv_get_unset();
  test_locv_gi();
  test_locv_ti();
  test_locv_invalid_typ();
  test_locv_negative_index();
  test_locv_high_index();
  test_locv_null_context();
  test_locv_set_get_roundtrip();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
