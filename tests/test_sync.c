#include "mythread/mythread.h"
#include <stdio.h>
#include <string.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

/*
 * These tests set up the global state so that Wait functions return
 * immediately (the waited-on state already matches the expected value).
 * The timeout exit() path in the Wait macros is avoided by setting up
 * correct state before calling Wait.
 */

/* TC-S01: sSetState / sWaitState roundtrip */
static void test_sm_set_wait(void) {
  TEST("sync-sm-set-wait");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  int saved_md_state = md.state;

  /* sSetState sets ti->state */
  sSetState(42);
  CHK(ti->state == 42, "sSetState did not set ti->state");

  /* sWaitState waits for md.state >= val. Set md.state to satisfy. */
  md.state = 100;
  sWaitState(50);  /* 100 >= 50, returns immediately */
  CHK(ti->state == 42, "sWaitState should not modify ti->state");

  ti = saved_ti;
  md.state = saved_md_state;
  PASS();
}

/* TC-S02: sWaitStater returns immediately when condition met */
static void test_sm_wait_stater(void) {
  TEST("sync-sm-wait-stater");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  int saved_md_state = md.state;
  md.state = 10;
  sWaitStater(20);  /* 10 <= 20, returns immediately */

  ti = saved_ti;
  md.state = saved_md_state;
  PASS();
}

/* TC-S03: mSetSubs / mWaitSubs roundtrip */
static void test_ms_set_wait(void) {
  TEST("sync-ms-set-wait");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  int saved_md_state = md.state;
  short saved_Nthreads = md.Nthreads;
  THREADINFO saved_threads[MCOREPC];
  memcpy(saved_threads, md.threads, sizeof(saved_threads));

  /* mSetSubs sets md.state */
  mSetSubs(5);
  CHK(md.state == 5, "mSetSubs did not set md.state");

  /* mWaitSubs checks all md.threads[i].state >= val.
   * Set Nthreads=1 so the loop iterates 0 times (no workers to check). */
  md.Nthreads = 1;
  mWaitSubs(10);  /* no workers, returns immediately */

  memcpy(md.threads, saved_threads, sizeof(saved_threads));
  md.Nthreads = saved_Nthreads;
  md.state = saved_md_state;
  ti = saved_ti;
  PASS();
}

/* TC-S04: mWaitSubsr returns immediately */
static void test_ms_wait_subsr(void) {
  TEST("sync-ms-wait-subsr");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  short saved_Nthreads = md.Nthreads;
  md.Nthreads = 1;  /* no workers to check */
  mWaitSubsr(10);

  md.Nthreads = saved_Nthreads;
  ti = saved_ti;
  PASS();
}

/* TC-S05: SG — gSetMain / gWaitMain roundtrip */
static void test_gm_set_wait(void) {
  TEST("sync-gm-set-wait");
  THREADINFO *saved_ti = ti;
  threadGroup *saved_gi = gi;
  THREADINFO t;
  threadGroup g;
  memset(&t, 0, sizeof(t));
  memset(&g, 0, sizeof(g));
  t.pg = &g;
  t.igrp = 0;
  ti = &t;
  gi = &g;

  /* gSetMain sets md.gstate[igrp*MSB] */
  int saved_gstate = md.gstate[0];
  gSetMain(7);
  CHK(md.gstate[0] == 7, "gSetMain did not set gstate");

  /* gWaitMain waits for gi->gstate >= val. Set to satisfy. */
  gi->gstate = 10;
  gWaitMain(5);  /* 10 >= 5, returns immediately */

  md.gstate[0] = saved_gstate;
  gi = saved_gi;
  ti = saved_ti;
  PASS();
}

/* TC-S06: MG — mSetGrps / mWaitGrps roundtrip */
static void test_mg_set_wait(void) {
  TEST("sync-mg-set-wait");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  int saved_md_state = md.state;
  short saved_ngrp = md.ngrp;
  threadGroup *saved_grps[MCLUST];
  memcpy(saved_grps, md.grps, sizeof(saved_grps));

  /* mSetGrps sets md.state and all group gstates */
  md.ngrp = 0;  /* no groups to iterate */
  mSetGrps(3);
  CHK(md.state == 3, "mSetGrps did not set md.state");

  /* mWaitGrps checks all group gstates */
  mWaitGrps(5);  /* no groups, returns immediately */

  memcpy(md.grps, saved_grps, sizeof(saved_grps));
  md.ngrp = saved_ngrp;
  md.state = saved_md_state;
  ti = saved_ti;
  PASS();
}

/* TC-S07: SG — sSetGrp / sWaitGrp roundtrip */
static void test_sg_set_wait(void) {
  TEST("sync-sg-set-wait");
  THREADINFO *saved_ti = ti;
  threadGroup *saved_gi = gi;
  THREADINFO t;
  threadGroup g;
  memset(&t, 0, sizeof(t));
  memset(&g, 0, sizeof(g));
  t.pg = &g;
  t.ind = 0;
  t.sib = 0;
  t.sie = 0;
  t.Nthreads = 1;
  ti = &t;
  gi = &g;

  gi->state = 100;
  sWaitGrp(50);  /* 100 >= 50, returns immediately */

  sSetGrp(99);   /* sets sstate[0] to 99, then waits for sib..sie (none when sie==sib) */

  gi = saved_gi;
  ti = saved_ti;
  PASS();
}

/* TC-S08: sWaitGrpr returns immediately when condition met */
static void test_sg_wait_grpr(void) {
  TEST("sync-sg-wait-grpr");
  THREADINFO *saved_ti = ti;
  threadGroup *saved_gi = gi;
  THREADINFO t;
  threadGroup g;
  memset(&t, 0, sizeof(t));
  memset(&g, 0, sizeof(g));
  t.pg = &g;
  ti = &t;
  gi = &g;

  gi->state = 5;
  sWaitGrpr(10);  /* 5 <= 10, returns immediately */

  gi = saved_gi;
  ti = saved_ti;
  PASS();
}

/* TC-S09: gWaitMainr returns immediately when condition met */
static void test_gm_wait_mainr(void) {
  TEST("sync-gm-wait-mainr");
  THREADINFO *saved_ti = ti;
  threadGroup *saved_gi = gi;
  THREADINFO t;
  threadGroup g;
  memset(&t, 0, sizeof(t));
  memset(&g, 0, sizeof(g));
  t.pg = &g;
  ti = &t;
  gi = &g;

  gi->gstate = 5;
  gWaitMainr(10);  /* 5 <= 10, returns immediately */

  gi = saved_gi;
  ti = saved_ti;
  PASS();
}

/* TC-S10: mWaitGrpsr returns immediately */
static void test_mg_wait_grpsr(void) {
  TEST("sync-mg-wait-grpsr");
  THREADINFO *saved_ti = ti;
  THREADINFO t;
  memset(&t, 0, sizeof(t));
  ti = &t;

  short saved_ngrp = md.ngrp;
  md.ngrp = 0;  /* no groups to check */
  mWaitGrpsr(10);

  md.ngrp = saved_ngrp;
  ti = saved_ti;
  PASS();
}

/* TC-S11: GetVInt returns the value pointed to */
static void test_getvint(void) {
  TEST("sync-getvint");
  volatile int v = 123;
  volatile int *vp = &v;
  CHK(GetVInt(vp) == 123, "GetVInt returned wrong value");

  v = -5;
  CHK(GetVInt(vp) == -5, "GetVInt returned wrong value after change");
  PASS();
}

/* TC-S12: GS — gSetSubs / gWaitSubs (with single worker) */
static void test_gs_set_wait(void) {
  TEST("sync-gs-set-wait");
  THREADINFO *saved_ti = ti;
  threadGroup *saved_gi = gi;
  THREADINFO t;
  threadGroup g;
  memset(&t, 0, sizeof(t));
  memset(&g, 0, sizeof(g));
  t.pg = &g;
  t.ind = 0;
  t.sib = 0;
  t.sie = 0;
  t.Nthreads = 1;
  ti = &t;
  gi = &g;

  gSetSubs(5);
  CHK(gi->state == 5, "gSetSubs did not set gi->state");

  /* Nthreads=1: sib=0, sie=0 → no workers to check in gWaitSubs */
  gWaitSubs(10);  /* returns immediately */

  gi = saved_gi;
  ti = saved_ti;
  PASS();
}

int main(void) {
  printf("\n=== mythread sync tests ===\n\n");

  test_sm_set_wait();
  test_sm_wait_stater();
  test_ms_set_wait();
  test_ms_wait_subsr();
  test_gm_set_wait();
  test_mg_set_wait();
  test_sg_set_wait();
  test_sg_wait_grpr();
  test_gm_wait_mainr();
  test_mg_wait_grpsr();
  test_getvint();
  test_gs_set_wait();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
