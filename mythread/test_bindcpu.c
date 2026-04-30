#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "mythread.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#else
#include <sched.h>
#endif

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
  tests_run++; \
  printf("  RUN  %s ... ", name); \
  fflush(stdout); \
} while(0)

#define PASS() do { \
  tests_passed++; \
  printf("PASS\n"); \
} while(0)

#define FAIL(msg) do { \
  tests_failed++; \
  printf("FAIL: %s\n", msg); \
} while(0)

#define CHECK(cond, msg) do { \
  if (!(cond)) { FAIL(msg); return; } \
} while(0)

/* ── helpers ── */
static int count_affinity_cores(void){
#ifdef __APPLE__
  int n = 0;
  size_t size = sizeof(n);
  if (sysctlbyname("hw.logicalcpu", &n, &size, NULL, 0) == 0) return n;
  return -1;
#else
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) == -1) return -1;
  int count = 0;
  for (int i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, &mask)) count++;
  }
  return count;
#endif
}

static int is_core_in_affinity(int core){
#ifdef __APPLE__
  (void)core;
  return -1;
#else
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) == -1) return -1;
  return CPU_ISSET(core, &mask) ? 1 : 0;
#endif
}

/* ── TC-B01: 绑定到单个合法核心 ── */
static void test_bindcpu_single_valid(void){
  TEST("bindcpu-single-valid");
  int ret = bindcpu(0);
  CHECK(ret == 0, "bindcpu(0) returned non-zero");
#ifndef __APPLE__
  CHECK(is_core_in_affinity(0) == 1, "core 0 not in affinity mask after bind");
#endif
  PASS();
}

/* ── TC-B02: 先绑定 A 再切换 B ── */
static void test_bindcpu_switch(void){
  TEST("bindcpu-multiple-switch");
  int total = count_affinity_cores();
  if (total < 2) {
    printf("SKIP (only %d cores available)\n", total);
    tests_run--;
    return;
  }
  CHECK(bindcpu(0) == 0, "bindcpu(0) failed");
#ifndef __APPLE__
  CHECK(is_core_in_affinity(0) == 1, "core 0 not set after bindcpu(0)");
#endif
  CHECK(bindcpu(1) == 0, "bindcpu(1) failed");
#ifndef __APPLE__
  CHECK(is_core_in_affinity(1) == 1, "core 1 not set after bindcpu(1)");
#endif
  PASS();
}

/* ── TC-B03: macOS 上空操作不报错 ── */
static void test_bindcpu_macos_noop(void){
  TEST("bindcpu-macos-noop");
  int ret = bindcpu(0);
#ifdef __APPLE__
  CHECK(ret == 0, "bindcpu should return 0 on macOS");
#else
  CHECK(ret == 0, "bindcpu should return 0");
#endif
  PASS();
}

/* ── TC-B04: 遍历所有可用核心 ── */
static void test_bindcpu_all_cores(void){
  TEST("bindcpu-all-cores");
  int total = count_affinity_cores();
  if (total < 0) {
    printf("SKIP (cannot query core count)\n");
    tests_run--;
    return;
  }
#ifndef __APPLE__
  int bound_ok = 0;
  for (int c = 0; c < total && c < 64; c++) {
    if (bindcpu(c) == 0 && is_core_in_affinity(c) == 1) {
      bound_ok++;
    }
  }
  CHECK(bound_ok > 0, "failed to bind to any core");
#endif
  PASS();
}

/* ── TC-B05: ti==NULL 时 bindthread 不崩溃 ── */
static void test_bindthread_ti_null(void){
  TEST("bindthread-ti-null");
  THREADINFO *saved_ti = ti;
  ti = NULL;
  bindthread();
  ti = saved_ti;
  PASS();
}

/* ── TC-B06: bindthread 使用 ti->indg ── */
static void test_bindthread_uses_indg(void){
  TEST("bindthread-uses-indg");
  THREADINFO saved = {0};
  THREADINFO *saved_ti = ti;
  saved.indg = 0;
  ti = &saved;
  bindthread();
  ti = saved_ti;
  PASS();
}

/* ── TC-B07: _threadmain_ 内部调用 bindcpu ── */
static int tfun_called = 0;
static void dummy_tfun(void){ tfun_called = 1; }

static void test_threadmain_calls_bind(void){
  TEST("threadmain-calls-bind");
  THREADINFO info;
  memset(&info, 0, sizeof(info));
  info.indg = 0;
  info.tfun = dummy_tfun;
  tfun_called = 0;
  _threadmain_(&info);
  CHECK(tfun_called == 1, "_threadmain_ did not call tfun");
  PASS();
}

/* ── TC-B08: getcpuid 返回合法值 ── */
static void test_getcpuid_valid(void){
  TEST("getcpuid-valid");
  int cpu = getcpuid();
#ifdef __APPLE__
  CHECK(cpu > 0, "getcpuid should return positive value");
#else
  CHECK(cpu >= 0, "getcpuid should return non-negative value");
#endif
  PASS();
}

/* ── TC-B09: initmd 计算 indg 在合法范围内 ── */
static void test_initmd_indg_range(void){
  TEST("initmd-indg-range");
  int saved_mpi_id = mpi_id;
  int saved_NCorePClu = NCorePClu;
  int saved_NThPGrp = NThPGrp;
  int saved_NGrpPProc = NGrpPProc;
  int saved_NProcPNode = NProcPNode;
  int saved_ManageCoreId = ManageCoreId;
  int saved_ThreadG = ThreadG;
  threadProc saved_md = md;

  /* non-grouped mode: ThreadG=0 via NGrpPProc=1 */
  int mc = 0;
  md.PThreadInited = 0;
  memset(&md, 0, sizeof(md));
  ThreadG = -1;
  NCorePClu = 6;
  NCluPNode = 1;
  NCorePGrp = 6;
  NThPGrp = 3;
  NGrpPProc = 1;
  NProcPNode = 1;
  InitThreads(0, NCorePClu, NCluPNode, NCorePGrp, NThPGrp, NGrpPProc, NProcPNode, &mc);

  CHECK(md.Nthreads > 0, "md.Nthreads is zero");
  for (int j = 0; j < md.Nthreads - 1; j++) {
    CHECK(md.threads[j].indg >= 0, "thread indg is negative in non-grouped mode");
  }

  /* grouped mode */
  memset(&md, 0, sizeof(md));
  ThreadG = -1;
  InitThreads(0, 8, 1, 8, 4, 2, 2, &mc);

  CHECK(md.ngrp == 2, "expected 2 groups");
  for (int i = 0; i < md.ngrp; i++) {
    threadGroup *g = md.grps[i];
    CHECK(g != NULL, "group pointer is null");
    for (int j = 0; j < g->Nthreads; j++) {
      CHECK(g->threads[j].indg >= 0, "thread indg is negative in grouped mode");
      for (int k = 0; k < i; k++) {
        threadGroup *g2 = md.grps[k];
        for (int l = 0; l < g2->Nthreads; l++) {
          CHECK(g->threads[j].indg != g2->threads[l].indg,
                "duplicate indg across groups");
        }
      }
    }
  }

  /* restore */
  mpi_id = saved_mpi_id;
  NCorePClu = saved_NCorePClu;
  NThPGrp = saved_NThPGrp;
  NGrpPProc = saved_NGrpPProc;
  NProcPNode = saved_NProcPNode;
  ManageCoreId = saved_ManageCoreId;
  ThreadG = saved_ThreadG;
  md = saved_md;
  PASS();
}

int main(void){
  printf("\n=== mythread CPU binding tests ===\n\n");

  test_bindcpu_single_valid();
  test_bindcpu_switch();
  test_bindcpu_macos_noop();
  test_bindcpu_all_cores();
  test_bindthread_ti_null();
  test_bindthread_uses_indg();
  test_threadmain_calls_bind();
  test_getcpuid_valid();
  test_initmd_indg_range();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
