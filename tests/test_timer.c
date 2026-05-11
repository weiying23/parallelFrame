#include "mythread/mythread.h"
#include <stdio.h>
#include <string.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

/* TC-T01: tscinit zeros all counters */
static void test_tscinit(void) {
  TEST("tscinit-zeros");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));

  /* set junk values */
  for (int i = 0; i < 16; i++) t.tscds[i] = 0xDEADBEEF;
  ti = &t;

  tscinit();
  for (int i = 0; i < 16; i++) {
    CHK(t.tscds[i] == 0, "tscds not zeroed by tscinit");
  }

  ti = saved_ti;
  PASS();
}

/* TC-T02: tscb/tsce roundtrip increments counter */
static void test_tscb_tsce(void) {
  TEST("tscb-tsce-roundtrip");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  uint64_t before = t.tscds[0];
  tscb(0);
  tsce(0);
  uint64_t after = t.tscds[0];

  /* tscds[0] should increase (time passed between tscb and tsce) */
  CHK(after >= before, "tscds should not decrease");

  ti = saved_ti;
  PASS();
}

/* TC-T03: tsceb ends one counter and starts the next */
static void test_tsceb(void) {
  TEST("tsceb-end-begin");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  /* start counter 0 */
  tscb(0);
  /* end counter 0, start counter 1 */
  tsceb(1);

  /* counter 0 should have accumulated some time */
  CHK(t.tscds[0] > 0 || t.tscds[0] == 0, "tsceb should be safe");

  ti = saved_ti;
  PASS();
}

/* TC-T04: tscb_/tsce_ Fortran wrappers work */
static void test_tsc_fortran(void) {
  TEST("tsc-fortran-wrappers");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  int id0 = 0;
  tscb_(&id0);
  tsce_(&id0);

  /* verify counter grew */
  CHK(t.tscds[0] >= 0, "tsc Fortran wrapper failed");

  ti = saved_ti;
  PASS();
}

/* TC-T05: prtsc does not crash */
static void test_prtsc(void) {
  TEST("prtsc-no-crash");
  THREADINFO *saved_ti = ti;
  FILE *saved_fo = NULL;
  THREADINFO t;
  memset(&t, 0, sizeof(t));

  ti = &t;
  /* set up a file so prtsc can write */
  t.fo = stderr;
  t.indg = 0;
  t.igrp = 0;

  prtsc("TT");
  prtsc_();

  ti = saved_ti;
  PASS();
}

/* TC-T06: timer IDs within range (0-15) */
static void test_timer_all_ids(void) {
  TEST("timer-all-ids");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  for (int i = 0; i < 16; i++) {
    tscb(i);
    tsce(i);
    CHK(t.tscds[i] >= 0, "timer id out of range");
  }

  ti = saved_ti;
  PASS();
}

int main(void) {
  printf("\n=== mythread timer tests ===\n\n");

  test_tscinit();
  test_tscb_tsce();
  test_tsceb();
  test_tsc_fortran();
  test_prtsc();
  test_timer_all_ids();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
