#include "mythread.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

/* ── TC-F01: Y_ONLY 分解 + 分配 + 释放 ── */
static void test_alloc_free_y_only(void) {
  TEST("alloc-free-y-only");
  mythread_decomp *dc = mythread_decomp_create(100, 200, 1, 4, 3,
      MYTHREAD_DECOMP_Y_ONLY);
  CHK(dc != NULL, "decomp_create failed");

  GroupField *gfs = calloc((size_t)dc->n_groups, sizeof(GroupField));
  CHK(gfs != NULL, "calloc gfs failed");

  for (int g = 0; g < dc->n_groups; g++) {
    int rc = group_alloc_field(&gfs[g], g, dc);
    CHK(rc == 0, "group_alloc_field failed");
    CHK(gfs[g].u_prev != NULL, "u_prev NULL");
    CHK(gfs[g].u_curr != NULL, "u_curr NULL");
    CHK(gfs[g].u_next != NULL, "u_next NULL");
    CHK(gfs[g].stride == dc->group_tiles[g].nx + 2 * dc->halo, "stride mismatch");
  }

  group_field_link_buffers(gfs, dc);

  /* 验证 Y 方向 send_to 指针互连 */
  for (int g = 0; g < dc->n_groups; g++) {
    int down = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_DOWN);
    if (down >= 0) {
      CHK(gfs[g].send_to_down == gfs[down].recv_from_up,
          "send_to_down not linked to neighbor recv_from_up");
    }
    int up = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_UP);
    if (up >= 0) {
      CHK(gfs[g].send_to_up == gfs[up].recv_from_down,
          "send_to_up not linked to neighbor recv_from_down");
    }
  }

  for (int g = 0; g < dc->n_groups; g++) group_free_field(&gfs[g]);
  free(gfs);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F02: XY_2D 分解 + 分配 + X 方向缓冲 ── */
static void test_alloc_free_xy2d(void) {
  TEST("alloc-free-xy2d");
  mythread_decomp *dc = mythread_decomp_create(160, 120, 1, 4, 2,
      MYTHREAD_DECOMP_XY_2D);
  CHK(dc != NULL, "decomp_create failed");

  GroupField *gfs = calloc((size_t)dc->n_groups, sizeof(GroupField));
  CHK(gfs != NULL, "calloc gfs failed");

  for (int g = 0; g < dc->n_groups; g++) {
    CHK(group_alloc_field(&gfs[g], g, dc) == 0, "group_alloc_field failed");
  }

  group_field_link_buffers(gfs, dc);

  /* 验证 X 方向缓冲存在（内部组应有 left/right 邻居） */
  int has_left = 0, has_right = 0;
  for (int g = 0; g < dc->n_groups; g++) {
    if (mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_LEFT) >= 0) {
      CHK(gfs[g].recv_from_left != NULL, "XY_2D recv_from_left NULL");
      has_left++;
    }
    if (mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_RIGHT) >= 0) {
      CHK(gfs[g].recv_from_right != NULL, "XY_2D recv_from_right NULL");
      has_right++;
    }
  }
  CHK(has_left > 0, "no group has left neighbor in XY_2D");
  CHK(has_right > 0, "no group has right neighbor in XY_2D");

  for (int g = 0; g < dc->n_groups; g++) group_free_field(&gfs[g]);
  free(gfs);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F03: GFIDX 宏线性索引 ── */
static void test_gfidx(void) {
  TEST("gfidx-linear-index");
  mythread_decomp *dc = mythread_decomp_create(10, 8, 2, 1, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  GroupField gf;
  CHK(group_alloc_field(&gf, 0, dc) == 0, "alloc failed");

  /* 写一个已知值验证索引 */
  gf.u_curr[GFIDX(&gf, 3, 5)] = 42.0;
  CHK(gf.u_curr[GFIDX(&gf, 3, 5)] == 42.0, "GFIDX readback mismatch");

  /* 二维排布：(y) * stride + x */
  size_t expected = (size_t)3 * (size_t)gf.stride + (size_t)5;
  CHK(GFIDX(&gf, 3, 5) == expected, "GFIDX formula mismatch");

  group_free_field(&gf);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F04: field_fill + field_swap ── */
static void test_fill_swap(void) {
  TEST("fill-and-swap");
  mythread_decomp *dc = mythread_decomp_create(10, 10, 1, 1, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  GroupField gf;
  group_alloc_field(&gf, 0, dc);

  group_field_fill(&gf, 1.5);
  CHK(gf.u_prev[GFIDX(&gf, 1, 1)] == 1.5, "fill prev failed");
  CHK(gf.u_curr[GFIDX(&gf, 2, 3)] == 1.5, "fill curr failed");

  /* 修改 curr，swap 后值应出现在 prev */
  gf.u_curr[GFIDX(&gf, 0, 0)] = 99.0;
  group_field_swap(&gf);
  CHK(gf.u_prev[GFIDX(&gf, 0, 0)] == 99.0,  "swap: old curr not in prev");
  CHK(gf.u_curr[GFIDX(&gf, 0, 0)] == 1.5,   "swap: old next not in curr");

  group_free_field(&gf);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F05: NUMA 节点信息输出 ── */
static void test_numa_info(void) {
  TEST("numa-info");
  mythread_decomp *dc = mythread_decomp_create(20, 20, 1, 2, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  GroupField gfs[2];
  for (int g = 0; g < 2; g++) {
    CHK(group_alloc_field(&gfs[g], g, dc) == 0, "alloc failed");
    /* numa_node 应为 -1（无 libnuma/macOS）或 ≥0（Linux+numa） */
    CHK(gfs[g].numa_node >= -1, "numa_node out of range");
  }
  group_field_dump(&gfs[0], 0, stdout);
  for (int g = 0; g < 2; g++) group_free_field(&gfs[g]);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F06: 重复 alloc/free 不泄漏 ── */
static void test_repeated_alloc_free(void) {
  TEST("repeated-alloc-free");
  mythread_decomp *dc = mythread_decomp_create(30, 30, 1, 2, 1,
      MYTHREAD_DECOMP_Y_ONLY);
  for (int iter = 0; iter < 100; iter++) {
    GroupField gf;
    CHK(group_alloc_field(&gf, iter % 2, dc) == 0, "alloc failed in loop");
    group_free_field(&gf);
  }
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F07: 域边界组分配 MPI halo 缓冲 ── */
static void test_boundary_mpi_buffers(void) {
  TEST("boundary-mpi-buffers");
  mythread_decomp *dc = mythread_decomp_create(50, 100, 2, 4, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  GroupField *gfs = calloc((size_t)dc->n_groups, sizeof(GroupField));

  for (int g = 0; g < dc->n_groups; g++)
    group_alloc_field(&gfs[g], g, dc);

  /* Y_ONLY: g=0 最下(y最小)→无DOWN邻居→需MPI下边界缓冲
             g=n-1 最上(y最大)→无UP邻居→需MPI上边界缓冲 */
  CHK(gfs[0].send_down != NULL && gfs[0].recv_down != NULL,
      "Y-lower boundary group (g=0) missing MPI down buffers");
  CHK(gfs[dc->n_groups - 1].send_up != NULL &&
      gfs[dc->n_groups - 1].recv_up != NULL,
      "Y-upper boundary group (g=n-1) missing MPI up buffers");

  /* 内部组不应有 Y 方向 MPI 缓冲 */
  for (int g = 1; g < dc->n_groups - 1; g++) {
    CHK(gfs[g].send_up == NULL && gfs[g].recv_up == NULL,
        "internal group has unexpected MPI up buffers");
    CHK(gfs[g].send_down == NULL && gfs[g].recv_down == NULL,
        "internal group has unexpected MPI down buffers");
  }

  for (int g = 0; g < dc->n_groups; g++) group_free_field(&gfs[g]);
  free(gfs);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F08: 单组模式（非分组兼容）── */
static void test_single_group(void) {
  TEST("single-group-mode");
  mythread_decomp *dc = mythread_decomp_create(60, 80, 1, 1, 3,
      MYTHREAD_DECOMP_Y_ONLY);
  GroupField gf;
  CHK(group_alloc_field(&gf, 0, dc) == 0, "alloc failed");
  CHK(gf.send_to_up == NULL && gf.send_to_down == NULL,
      "single group should have no intra-group neighbors");
  /* 单组 = 全域边界 → MPI 缓冲应全部分配 */
  CHK(gf.send_up && gf.send_down && gf.send_left && gf.send_right,
      "single group missing MPI buffers");
  group_free_field(&gf);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F09: intra halo 交换数据正确性 ── */
static void test_halo_intra_correctness(void) {
  TEST("halo-intra-correctness");
  /* 3 组 Y_ONLY，每组 4x3，halo=1 → 含 halo: 6x5 */
  mythread_decomp *dc = mythread_decomp_create(4, 9, 1, 3, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  GroupField *gfs = calloc((size_t)dc->n_groups, sizeof(GroupField));
  for (int g = 0; g < dc->n_groups; g++)
    CHK(group_alloc_field(&gfs[g], g, dc) == 0, "alloc failed");
  group_field_link_buffers(gfs, dc);

  /* 每组填充不同的值：gf[g] 全部填 g+1 */
  for (int g = 0; g < dc->n_groups; g++)
    group_field_fill(&gfs[g], (double)(g + 1));

  /* 执行组间 halo 交换 */
  mythread_halo_exchange_intra(gfs, dc);

  /* 验证：
     gf[0]: 全部=1.0，上halo(行0)应来自 gf[1]=2.0，下halo 无邻居=保持1.0
     gf[1]: 全部=2.0，上halo 来自 gf[2]=3.0，下halo 来自 gf[0]=1.0
     gf[2]: 全部=3.0，上halo 无邻居=保持3.0，下halo 来自 gf[1]=2.0
   */
  int halo = dc->halo;
  /* gf[0]: 上 halo 行应该变成 2.0（来自 gf[1] 的第一行内部数据=2.0） */
  CHK(gfs[0].u_curr[GFIDX(&gfs[0], 0, halo)] == 2.0,
      "gf[0] up-halo should be 2.0 from gf[1]");
  /* gf[0]: 下 halo 行无人填充，保持 1.0 */
  CHK(gfs[0].u_curr[GFIDX(&gfs[0], dc->group_tiles[0].ny + halo, halo)] == 1.0,
      "gf[0] down-halo should stay 1.0 (no neighbor)");

  /* gf[1]: 下 halo 行来自 gf[0] 的数据=1.0 */
  CHK(gfs[1].u_curr[GFIDX(&gfs[1], dc->group_tiles[1].ny + halo, halo)] == 1.0,
      "gf[1] down-halo should be 1.0 from gf[0]");
  /* gf[1]: 上 halo 行来自 gf[2] 的数据=3.0 */
  CHK(gfs[1].u_curr[GFIDX(&gfs[1], 0, halo)] == 3.0,
      "gf[1] up-halo should be 3.0 from gf[2]");

  /* gf[2]: 下 halo 来自 gf[1]=2.0 */
  CHK(gfs[2].u_curr[GFIDX(&gfs[2], dc->group_tiles[2].ny + halo, halo)] == 2.0,
      "gf[2] down-halo should be 2.0 from gf[1]");

  for (int g = 0; g < dc->n_groups; g++) group_free_field(&gfs[g]);
  free(gfs);
  mythread_decomp_free(dc);
  PASS();
}

/* ── TC-F10: intra halo XY_2D 四方向验证 ── */
static void test_halo_intra_xy2d(void) {
  TEST("halo-intra-xy2d");
  /* 4 组 2x2，每组 3x3 内部，halo=1 → 含 halo: 5x5 */
  mythread_decomp *dc = mythread_decomp_create(6, 6, 1, 4, 0,
      MYTHREAD_DECOMP_XY_2D);
  GroupField *gfs = calloc((size_t)dc->n_groups, sizeof(GroupField));
  for (int g = 0; g < dc->n_groups; g++)
    CHK(group_alloc_field(&gfs[g], g, dc) == 0, "alloc failed");
  group_field_link_buffers(gfs, dc);

  /* 每组填不同值 */
  for (int g = 0; g < dc->n_groups; g++)
    group_field_fill(&gfs[g], (double)(g + 10));

  mythread_halo_exchange_intra(gfs, dc);

  int halo = dc->halo;
  /* XY_2D 2x2 grid: g0(左下), g1(右下), g2(左上), g3(右上)
     g0: right=1(右), up=2(上), down=-1, left=-1
       → 右 halo 来自 g1=11.0, 上 halo(row 0) 来自 g2=12.0 */
  CHK(gfs[0].u_curr[GFIDX(&gfs[0], halo + 1, dc->group_tiles[0].nx + halo)] == 11.0,
      "XY_2D gf[0] right-halo should be 11.0 from gf[1]");
  CHK(gfs[0].u_curr[GFIDX(&gfs[0], 0, halo)] == 12.0,
      "XY_2D gf[0] up-halo(row 0) should be 12.0 from gf[2]");

  /* g3(右上): left=2(左), down=1(下), up=-1, right=-1
       → 左 halo 来自 g2=12.0, 下 halo(row ny+halo) 来自 g1=11.0 */
  CHK(gfs[3].u_curr[GFIDX(&gfs[3], halo, 0)] == 12.0,
      "XY_2D gf[3] left-halo should be 12.0 from gf[2]");
  CHK(gfs[3].u_curr[GFIDX(&gfs[3], dc->group_tiles[3].ny + halo, halo)] == 11.0,
      "XY_2D gf[3] down-halo should be 11.0 from gf[1]");

  for (int g = 0; g < dc->n_groups; g++) group_free_field(&gfs[g]);
  free(gfs);
  mythread_decomp_free(dc);
  PASS();
}

int main(void) {
  printf("\n=== mythread GroupField + Halo tests ===\n\n");

  test_alloc_free_y_only();
  test_alloc_free_xy2d();
  test_gfidx();
  test_fill_swap();
  test_numa_info();
  test_repeated_alloc_free();
  test_boundary_mpi_buffers();
  test_single_group();
  test_halo_intra_correctness();
  test_halo_intra_xy2d();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
