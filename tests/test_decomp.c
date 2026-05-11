#include "mythread/mythread_decomp.h"
#include <stdio.h>
#include <string.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

/* TC-D01: create/free Y_ONLY basic */
static void test_decomp_y_only_basic(void) {
  TEST("decomp-y-only-basic");
  mythread_decomp *dc = mythread_decomp_create(100, 80, 1, 4, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  CHK(dc != NULL, "create failed");
  CHK(dc->n_groups == 4, "n_groups mismatch");
  CHK(dc->gx == 1, "gx should be 1 for Y_ONLY");
  CHK(dc->gy == 4, "gy should be n_groups for Y_ONLY");
  CHK(dc->group_tiles != NULL, "group_tiles NULL");

  /* all groups have same nx=100 */
  for (int g = 0; g < 4; g++) {
    CHK(dc->group_tiles[g].nx == 100, "nx mismatch in Y_ONLY");
    CHK(dc->group_tiles[g].ny > 0, "ny is zero");
    CHK(dc->group_tiles[g].x_begin == 0, "x_begin should be 0");
    CHK(dc->group_tiles[g].x_end == 100, "x_end should be 100");
  }

  /* verify y ranges partition [0, 80) exactly */
  int y_sum = 0;
  for (int g = 0; g < 4; g++) y_sum += dc->group_tiles[g].ny;
  CHK(y_sum == 80, "y partition does not sum to domain total");

  /* verify contiguity */
  for (int g = 1; g < 4; g++) {
    CHK(dc->group_tiles[g].y_begin == dc->group_tiles[g-1].y_end,
        "y partition not contiguous");
  }

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D02: Y_ONLY neighbor topology */
static void test_decomp_y_only_neighbors(void) {
  TEST("decomp-y-only-neighbors");
  mythread_decomp *dc = mythread_decomp_create(10, 30, 1, 3, 0,
      MYTHREAD_DECOMP_Y_ONLY);

  /* g=0: no down, up=1 */
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_DOWN) == -1,
      "g0 should have no down");
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_UP) == 1,
      "g0 up should be 1");
  CHK(mythread_decomp_is_domain_boundary(dc, 0, MYTHREAD_NEIGHBOR_DOWN),
      "g0 down should be domain boundary");
  CHK(!mythread_decomp_is_domain_boundary(dc, 0, MYTHREAD_NEIGHBOR_UP),
      "g0 up should not be domain boundary");

  /* g=1: down=0, up=2 */
  CHK(mythread_decomp_neighbor(dc, 1, MYTHREAD_NEIGHBOR_DOWN) == 0,
      "g1 down should be 0");
  CHK(mythread_decomp_neighbor(dc, 1, MYTHREAD_NEIGHBOR_UP) == 2,
      "g1 up should be 2");

  /* g=2: down=1, no up */
  CHK(mythread_decomp_neighbor(dc, 2, MYTHREAD_NEIGHBOR_DOWN) == 1,
      "g2 down should be 1");
  CHK(mythread_decomp_neighbor(dc, 2, MYTHREAD_NEIGHBOR_UP) == -1,
      "g2 should have no up");
  CHK(mythread_decomp_is_domain_boundary(dc, 2, MYTHREAD_NEIGHBOR_UP),
      "g2 up should be domain boundary");

  /* Y_ONLY: no left/right */
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_LEFT) == -1,
      "Y_ONLY should have no left");
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_RIGHT) == -1,
      "Y_ONLY should have no right");

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D03: XY_2D basic */
static void test_decomp_xy_2d_basic(void) {
  TEST("decomp-xy2d-basic");
  mythread_decomp *dc = mythread_decomp_create(120, 90, 1, 6, 0,
      MYTHREAD_DECOMP_XY_2D);
  CHK(dc != NULL, "create failed");
  CHK(dc->n_groups == 6, "n_groups mismatch");
  CHK(dc->gx * dc->gy == 6, "grid does not match n_groups");
  CHK(dc->gx > 1 || dc->gy > 1, "2D should have multi-dim grid");

  /* verify tile coverage */
  int area = 0;
  for (int g = 0; g < 6; g++)
    area += dc->group_tiles[g].nx * dc->group_tiles[g].ny;
  CHK(area == 120 * 90, "XY_2D area mismatch");

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D04: XY_2D neighbor topology */
static void test_decomp_xy_2d_neighbors(void) {
  TEST("decomp-xy2d-neighbors");
  mythread_decomp *dc = mythread_decomp_create(10, 10, 1, 4, 0,
      MYTHREAD_DECOMP_XY_2D);

  /* 4 groups: the grid will be 2x2 */
  int gx = dc->gx, gy = dc->gy;
  CHK(gx == 2 && gy == 2, "4 groups should map to 2x2 grid");

  /* g0: bottom-left */
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_DOWN) == -1,
      "g0 down");
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_LEFT) == -1,
      "g0 left");
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_UP) == 2,
      "g0 up");
  CHK(mythread_decomp_neighbor(dc, 0, MYTHREAD_NEIGHBOR_RIGHT) == 1,
      "g0 right");

  /* g1: bottom-right */
  CHK(mythread_decomp_neighbor(dc, 1, MYTHREAD_NEIGHBOR_DOWN) == -1,
      "g1 down");
  CHK(mythread_decomp_neighbor(dc, 1, MYTHREAD_NEIGHBOR_RIGHT) == -1,
      "g1 right");
  CHK(mythread_decomp_neighbor(dc, 1, MYTHREAD_NEIGHBOR_LEFT) == 0,
      "g1 left");
  CHK(mythread_decomp_neighbor(dc, 1, MYTHREAD_NEIGHBOR_UP) == 3,
      "g1 up");

  /* g2: top-left */
  CHK(mythread_decomp_neighbor(dc, 2, MYTHREAD_NEIGHBOR_UP) == -1,
      "g2 up");
  CHK(mythread_decomp_neighbor(dc, 2, MYTHREAD_NEIGHBOR_LEFT) == -1,
      "g2 left");
  CHK(mythread_decomp_neighbor(dc, 2, MYTHREAD_NEIGHBOR_DOWN) == 0,
      "g2 down");
  CHK(mythread_decomp_neighbor(dc, 2, MYTHREAD_NEIGHBOR_RIGHT) == 3,
      "g2 right");

  /* g3: top-right */
  CHK(mythread_decomp_neighbor(dc, 3, MYTHREAD_NEIGHBOR_UP) == -1,
      "g3 up");
  CHK(mythread_decomp_neighbor(dc, 3, MYTHREAD_NEIGHBOR_RIGHT) == -1,
      "g3 right");
  CHK(mythread_decomp_neighbor(dc, 3, MYTHREAD_NEIGHBOR_DOWN) == 1,
      "g3 down");
  CHK(mythread_decomp_neighbor(dc, 3, MYTHREAD_NEIGHBOR_LEFT) == 2,
      "g3 left");

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D05: worker tiles */
static void test_decomp_worker_tiles(void) {
  TEST("decomp-worker-tiles");
  mythread_decomp *dc = mythread_decomp_create(60, 40, 1, 2, 3,
      MYTHREAD_DECOMP_Y_ONLY);
  CHK(dc->n_workers_per_group == 3, "workers per group mismatch");

  /* each group gets 3 workers */
  for (int g = 0; g < 2; g++) {
    for (int w = 0; w < 3; w++) {
      const mythread_tile *wt = mythread_decomp_worker_tile(dc, g, w);
      CHK(wt != NULL, "worker tile NULL");
      CHK(wt->nx > 0 && wt->nx <= 60, "worker nx out of range");
      CHK(wt->ny > 0, "worker ny is zero");
    }
  }

  /* verify workers partition group tiles (area matches, since workers may split X or Y) */
  for (int g = 0; g < 2; g++) {
    const mythread_tile *gt = &dc->group_tiles[g];
    int area_sum = 0;
    for (int w = 0; w < 3; w++) {
      const mythread_tile *wt = mythread_decomp_worker_tile(dc, g, w);
      CHK(wt->y_begin >= gt->y_begin, "worker y_begin below group");
      CHK(wt->y_end <= gt->y_end, "worker y_end above group");
      CHK(wt->x_begin >= gt->x_begin, "worker x_begin left of group");
      CHK(wt->x_end <= gt->x_end, "worker x_end right of group");
      area_sum += wt->nx * wt->ny;
    }
    CHK(area_sum == gt->nx * gt->ny, "worker area partition does not sum to group area");
  }

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D06: invalid params return NULL */
static void test_decomp_invalid(void) {
  TEST("decomp-invalid-params");
  CHK(mythread_decomp_create(0, 10, 1, 1, 0, MYTHREAD_DECOMP_Y_ONLY) == NULL,
      "zero domain_nx accepted");
  CHK(mythread_decomp_create(10, 0, 1, 1, 0, MYTHREAD_DECOMP_Y_ONLY) == NULL,
      "zero domain_ny accepted");
  CHK(mythread_decomp_create(10, 10, -1, 1, 0, MYTHREAD_DECOMP_Y_ONLY) == NULL,
      "negative halo accepted");
  CHK(mythread_decomp_create(10, 10, 1, 0, 0, MYTHREAD_DECOMP_Y_ONLY) == NULL,
      "zero n_groups accepted");
  CHK(mythread_decomp_create(10, 10, 1, 1, -1, MYTHREAD_DECOMP_Y_ONLY) == NULL,
      "negative workers accepted");
  PASS();
}

/* TC-D07: decomp_free on NULL is safe */
static void test_decomp_free_null(void) {
  TEST("decomp-free-null");
  mythread_decomp_free(NULL);  /* should not crash */
  PASS();
}

/* TC-D08: decomp_dump does not crash */
static void test_decomp_dump(void) {
  TEST("decomp-dump");
  mythread_decomp *dc = mythread_decomp_create(20, 20, 1, 2, 2,
      MYTHREAD_DECOMP_XY_2D);
  /* dump with NULL out should not crash */
  mythread_decomp_dump(NULL, stdout);
  mythread_decomp_dump(dc, NULL);
  mythread_decomp_dump(dc, stdout);
  mythread_decomp_free(dc);
  PASS();
}

/* TC-D09: single group is domain boundary on all sides */
static void test_decomp_single_group_boundary(void) {
  TEST("decomp-single-group-boundary");
  mythread_decomp *dc = mythread_decomp_create(10, 10, 1, 1, 0,
      MYTHREAD_DECOMP_Y_ONLY);

  CHK(mythread_decomp_is_domain_boundary(dc, 0, MYTHREAD_NEIGHBOR_UP),
      "single group up should be boundary");
  CHK(mythread_decomp_is_domain_boundary(dc, 0, MYTHREAD_NEIGHBOR_DOWN),
      "single group down should be boundary");
  CHK(mythread_decomp_is_domain_boundary(dc, 0, MYTHREAD_NEIGHBOR_LEFT),
      "single group left should be boundary");
  CHK(mythread_decomp_is_domain_boundary(dc, 0, MYTHREAD_NEIGHBOR_RIGHT),
      "single group right should be boundary");

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D10: unequal domain partition handles remainder correctly */
static void test_decomp_uneven_partition(void) {
  TEST("decomp-uneven-partition");
  /* 11 / 3 = 3 + remainder 2 */
  mythread_decomp *dc = mythread_decomp_create(10, 11, 1, 3, 0,
      MYTHREAD_DECOMP_Y_ONLY);

  /* first rem groups get 1 extra */
  CHK(dc->group_tiles[0].ny == 3 + 1, "g0 ny");
  CHK(dc->group_tiles[1].ny == 3 + 1, "g1 ny");
  CHK(dc->group_tiles[2].ny == 3 + 0, "g2 ny");
  CHK(dc->group_tiles[0].ny + dc->group_tiles[1].ny + dc->group_tiles[2].ny == 11,
      "ny sum mismatch");

  mythread_decomp_free(dc);
  PASS();
}

/* TC-D11: nworkers_per_group=0 skips worker_tiles */
static void test_decomp_zero_workers(void) {
  TEST("decomp-zero-workers");
  mythread_decomp *dc = mythread_decomp_create(10, 10, 1, 4, 0,
      MYTHREAD_DECOMP_Y_ONLY);
  CHK(dc != NULL, "create failed");
  CHK(dc->worker_tiles == NULL, "worker_tiles should be NULL when nworkers=0");
  CHK(dc->n_workers_per_group == 0, "n_workers_per_group should be 0");
  mythread_decomp_free(dc);
  PASS();
}

int main(void) {
  printf("\n=== mythread decomp tests ===\n\n");

  test_decomp_y_only_basic();
  test_decomp_y_only_neighbors();
  test_decomp_xy_2d_basic();
  test_decomp_xy_2d_neighbors();
  test_decomp_worker_tiles();
  test_decomp_invalid();
  test_decomp_free_null();
  test_decomp_dump();
  test_decomp_single_group_boundary();
  test_decomp_uneven_partition();
  test_decomp_zero_workers();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
