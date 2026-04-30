#include "mythread_decomp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── 辅助：均匀切分区间 ── */
static void split_range(int begin, int end, int parts, int index,
                        int *sub_begin, int *sub_end) {
  int total = end - begin;
  int base = total / parts;
  int rem  = total % parts;
  int extra = (index < rem) ? 1 : 0;
  int offset = index * base + ((index < rem) ? index : rem);

  *sub_begin = begin + offset;
  *sub_end   = *sub_begin + base + extra;
}

/* ── 辅助：选择接近方形的 2D 分解 ── */
static void choose_2d_grid(int parts, int span_x, int span_y,
                           int *px, int *py) {
  int best_x = 1, best_y = 1;
  int best_penalty = 0;
  double best_balance = 0.0;
  int found = 0;

  if (parts <= 1) { *px = 1; *py = 1; return; }
  if (span_x < 1) span_x = 1;
  if (span_y < 1) span_y = 1;

  for (int x = 1; x <= parts; x++) {
    if (parts % x != 0) continue;
    int y = parts / x;
    int penalty = 0;
    if (x > span_x) penalty += x - span_x;
    if (y > span_y) penalty += y - span_y;
    double balance = fabs((double)span_x * y - (double)span_y * x);

    if (!found || penalty < best_penalty ||
        (penalty == best_penalty && balance < best_balance)) {
      found = 1;
      best_penalty = penalty;
      best_balance  = balance;
      best_x = x; best_y = y;
    }
  }
  *px = best_x; *py = best_y;
}

/* ── 构建组间邻居拓扑 ── */
static void build_group_neighbors(mythread_decomp *dc) {
  int gx = dc->gx, gy = dc->gy;

  for (int g = 0; g < dc->n_groups; g++) {
    int *nb = &dc->group_neighbors[g * MYTHREAD_NEIGHBOR_COUNT];

    if (dc->policy == MYTHREAD_DECOMP_Y_ONLY) {
      nb[MYTHREAD_NEIGHBOR_UP]    = (g > 0)                ? (g - 1) : -1;
      nb[MYTHREAD_NEIGHBOR_DOWN]  = (g < dc->n_groups - 1) ? (g + 1) : -1;
      nb[MYTHREAD_NEIGHBOR_LEFT]  = -1;
      nb[MYTHREAD_NEIGHBOR_RIGHT] = -1;
    } else {
      int gix = g % gx, giy = g / gx;
      nb[MYTHREAD_NEIGHBOR_UP]    = (giy > 0)      ? (g - gx) : -1;
      nb[MYTHREAD_NEIGHBOR_DOWN]  = (giy < gy - 1) ? (g + gx) : -1;
      nb[MYTHREAD_NEIGHBOR_LEFT]  = (gix > 0)      ? (g - 1)  : -1;
      nb[MYTHREAD_NEIGHBOR_RIGHT] = (gix < gx - 1) ? (g + 1)  : -1;
    }
  }
}

/* ── 公共 API ── */

mythread_decomp *mythread_decomp_create(int domain_nx, int domain_ny, int halo,
                                        int n_groups, int n_workers_per_group,
                                        mythread_decomp_policy policy) {
  if (domain_nx <= 0 || domain_ny <= 0 || halo < 0 ||
      n_groups <= 0 || n_workers_per_group < 0) {
    return NULL;
  }

  mythread_decomp *dc = calloc(1, sizeof(*dc));
  if (!dc) return NULL;

  dc->domain_nx           = domain_nx;
  dc->domain_ny           = domain_ny;
  dc->halo                = halo;
  dc->n_groups            = n_groups;
  dc->n_workers_per_group = n_workers_per_group;
  dc->policy              = policy;

  /* 组级 2D 网格 */
  if (policy == MYTHREAD_DECOMP_Y_ONLY) {
    dc->gx = 1;
    dc->gy = n_groups;
  } else {
    choose_2d_grid(n_groups, domain_nx, domain_ny, &dc->gx, &dc->gy);
  }

  /* 分配组级 tiles */
  dc->group_tiles = calloc((size_t)n_groups, sizeof(*dc->group_tiles));
  if (!dc->group_tiles) goto fail;

  for (int g = 0; g < n_groups; g++) {
    int gx_id = (policy == MYTHREAD_DECOMP_Y_ONLY) ? 0 : (g % dc->gx);
    int gy_id = (policy == MYTHREAD_DECOMP_Y_ONLY) ? g : (g / dc->gx);

    split_range(0, domain_nx, (policy == MYTHREAD_DECOMP_Y_ONLY) ? 1 : dc->gx,
                gx_id, &dc->group_tiles[g].x_begin, &dc->group_tiles[g].x_end);
    split_range(0, domain_ny, (policy == MYTHREAD_DECOMP_Y_ONLY) ? n_groups : dc->gy,
                gy_id, &dc->group_tiles[g].y_begin, &dc->group_tiles[g].y_end);

    dc->group_tiles[g].nx = dc->group_tiles[g].x_end - dc->group_tiles[g].x_begin;
    dc->group_tiles[g].ny = dc->group_tiles[g].y_end - dc->group_tiles[g].y_begin;
  }

  /* 组间邻居拓扑 */
  dc->group_neighbors = calloc((size_t)n_groups * MYTHREAD_NEIGHBOR_COUNT, sizeof(int));
  if (!dc->group_neighbors) goto fail;
  for (int i = 0; i < n_groups * MYTHREAD_NEIGHBOR_COUNT; i++)
    dc->group_neighbors[i] = -1;
  build_group_neighbors(dc);

  /* 组内 worker 分解 */
  if (n_workers_per_group > 0) {
    int total_workers = n_groups * n_workers_per_group;
    dc->worker_tiles = calloc((size_t)total_workers, sizeof(*dc->worker_tiles));
    if (!dc->worker_tiles) goto fail;

    for (int g = 0; g < n_groups; g++) {
      const mythread_tile *gt = &dc->group_tiles[g];

      /* worker 在 X+Y 方向均匀切分组的子域 */
      int wx = 1, wy = 1;
      choose_2d_grid(n_workers_per_group, gt->nx, gt->ny, &wx, &wy);

      for (int w = 0; w < n_workers_per_group; w++) {
        int widx = w % wx, widy = w / wx;
        mythread_tile *wt = &dc->worker_tiles[g * n_workers_per_group + w];

        split_range(gt->x_begin, gt->x_end, wx, widx, &wt->x_begin, &wt->x_end);
        split_range(gt->y_begin, gt->y_end, wy, widy, &wt->y_begin, &wt->y_end);
        wt->nx = wt->x_end - wt->x_begin;
        wt->ny = wt->y_end - wt->y_begin;
      }
    }
  }

  return dc;

fail:
  mythread_decomp_free(dc);
  return NULL;
}

void mythread_decomp_free(mythread_decomp *dc) {
  if (!dc) return;
  free(dc->group_tiles);
  free(dc->group_neighbors);
  free(dc->worker_tiles);
  memset(dc, 0, sizeof(*dc));
  free(dc);
}

void mythread_decomp_dump(const mythread_decomp *dc, FILE *out) {
  if (!dc || !out) return;

  fprintf(out, "=== domain decomposition ===\n");
  fprintf(out, "domain  : %d x %d, halo=%d\n", dc->domain_nx, dc->domain_ny, dc->halo);
  fprintf(out, "policy  : %s\n", dc->policy == MYTHREAD_DECOMP_Y_ONLY ? "Y_ONLY" : "XY_2D");
  fprintf(out, "groups  : %d (grid %d x %d), workers/group: %d\n",
          dc->n_groups, dc->gx, dc->gy, dc->n_workers_per_group);

  for (int g = 0; g < dc->n_groups; g++) {
    const mythread_tile *t = &dc->group_tiles[g];
    fprintf(out, "  group %d: x=[%d,%d) y=[%d,%d) nx=%d ny=%d  neighbors:",
            g, t->x_begin, t->x_end, t->y_begin, t->y_end, t->nx, t->ny);
    for (int d = 0; d < MYTHREAD_NEIGHBOR_COUNT; d++) {
      int nb = mythread_decomp_neighbor(dc, g, (mythread_neighbor_dir)d);
      if (nb >= 0) fprintf(out, " %s=%d", (const char*[]){"up","down","left","right"}[d], nb);
    }
    fprintf(out, "\n");

    for (int w = 0; w < dc->n_workers_per_group; w++) {
      const mythread_tile *wt = mythread_decomp_worker_tile(dc, g, w);
      fprintf(out, "    worker %d: x=[%d,%d) y=[%d,%d)\n",
              w, wt->x_begin, wt->x_end, wt->y_begin, wt->y_end);
    }
  }
  fprintf(out, "=== end decomposition ===\n");
}
