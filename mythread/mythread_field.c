#include "mythread_field.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── NUMA 支持检测 ── */
#if defined(__linux__) && !defined(NO_LIBNUMA)
# if __has_include(<numa.h>)
#  include <numa.h>
#  define HAS_LIBNUMA 1
# endif
#endif

#ifndef HAS_LIBNUMA
#define HAS_LIBNUMA 0
#endif

/* ── 平台无关的分配包装 ── */
static void *field_alloc(size_t bytes, int node) {
  if (bytes == 0) return NULL;
#if HAS_LIBNUMA
  if (node >= 0 && numa_available() >= 0) {
    return numa_alloc_onnode(bytes, node);
  }
#endif
  void *p = hmalloc(bytes);
  if (p) memset(p, 0, bytes);  /* 强制触发 page fault，确保 first-touch */
  return p;
}

static void field_free(void *p, size_t bytes) {
  if (!p) return;
  (void)bytes;
#if HAS_LIBNUMA
  if (numa_available() >= 0) {
    numa_free(p, bytes);
    return;
  }
#endif
  free(p);
}

/* ── NUMA 节点查询 ── */
static int get_numa_node(void) {
#if HAS_LIBNUMA
  if (numa_available() >= 0) {
    int cpu = sched_getcpu();
    return numa_node_of_cpu(cpu);
  }
#endif
  return -1;
}

/* ═══════════════════════════════════════════════════════════════ */

int group_alloc_field(GroupField *gf, int gid, const mythread_decomp *dc) {
  if (!gf || !dc || gid < 0 || gid >= dc->n_groups) return -1;

  const mythread_tile *tile = &dc->group_tiles[gid];
  int halo = dc->halo;
  int nx   = tile->nx;
  int ny   = tile->ny;
  int node = get_numa_node();

  memset(gf, 0, sizeof(*gf));

  gf->ny_padded   = ny + 2 * halo;
  gf->nx_padded   = nx + 2 * halo;
  gf->stride      = gf->nx_padded;
  gf->plane_bytes = (size_t)gf->ny_padded * (size_t)gf->nx_padded * sizeof(double);
  gf->numa_node   = node;
  gf->group_energy = 0.0;

  /* 波场平面 */
  gf->u_prev = (double*)field_alloc(gf->plane_bytes, node);
  gf->u_curr = (double*)field_alloc(gf->plane_bytes, node);
  gf->u_next = (double*)field_alloc(gf->plane_bytes, node);
  if (!gf->u_prev || !gf->u_curr || !gf->u_next) goto fail;

  /* 组间 Y 方向 halo 缓冲 */
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_UP) >= 0) {
    gf->recv_from_up = (double*)field_alloc((size_t)nx * sizeof(double), node);
    if (!gf->recv_from_up) goto fail;
  }
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_DOWN) >= 0) {
    gf->recv_from_down = (double*)field_alloc((size_t)nx * sizeof(double), node);
    if (!gf->recv_from_down) goto fail;
  }

  /* 组间 X 方向 halo 缓冲（仅 XY_2D 有邻居时分配） */
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_LEFT) >= 0) {
    gf->recv_from_left = (double*)field_alloc((size_t)ny * sizeof(double), node);
    if (!gf->recv_from_left) goto fail;
  }
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_RIGHT) >= 0) {
    gf->recv_from_right = (double*)field_alloc((size_t)ny * sizeof(double), node);
    if (!gf->recv_from_right) goto fail;
  }

  /* MPI halo 缓冲 — 仅域边界组 */
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_UP)) {
    gf->send_up = (double*)field_alloc((size_t)nx * sizeof(double), node);
    gf->recv_up = (double*)field_alloc((size_t)nx * sizeof(double), node);
    if (!gf->send_up || !gf->recv_up) goto fail;
  }
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_DOWN)) {
    gf->send_down = (double*)field_alloc((size_t)nx * sizeof(double), node);
    gf->recv_down = (double*)field_alloc((size_t)nx * sizeof(double), node);
    if (!gf->send_down || !gf->recv_down) goto fail;
  }
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_LEFT)) {
    gf->send_left  = (double*)field_alloc((size_t)ny * sizeof(double), node);
    gf->recv_left  = (double*)field_alloc((size_t)ny * sizeof(double), node);
    if (!gf->send_left || !gf->recv_left) goto fail;
  }
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_RIGHT)) {
    gf->send_right = (double*)field_alloc((size_t)ny * sizeof(double), node);
    gf->recv_right = (double*)field_alloc((size_t)ny * sizeof(double), node);
    if (!gf->send_right || !gf->recv_right) goto fail;
  }

  return 0;

fail:
  group_free_field(gf);
  return -1;
}

void group_free_field(GroupField *gf) {
  if (!gf) return;

  field_free(gf->u_prev,  gf->plane_bytes);
  field_free(gf->u_curr,  gf->plane_bytes);
  field_free(gf->u_next,  gf->plane_bytes);

  /* 组间 halo */
  field_free(gf->recv_from_up,    0);
  field_free(gf->recv_from_down,  0);
  field_free(gf->recv_from_left,  0);
  field_free(gf->recv_from_right, 0);
  /* send_to_* 不单独释放：它们是指向邻居 recv 的指针，由对应邻居释放 */

  /* MPI halo */
  field_free(gf->send_up,    0);
  field_free(gf->recv_up,    0);
  field_free(gf->send_down,  0);
  field_free(gf->recv_down,  0);
  field_free(gf->send_left,  0);
  field_free(gf->recv_left,  0);
  field_free(gf->send_right, 0);
  field_free(gf->recv_right, 0);

  memset(gf, 0, sizeof(*gf));
}

void group_field_link_buffers(GroupField *gfields, const mythread_decomp *dc) {
  if (!gfields || !dc) return;

  for (int g = 0; g < dc->n_groups; g++) {
    GroupField *gf = &gfields[g];

    int nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_UP);
    if (nb >= 0) {
      gf->send_to_up = gfields[nb].recv_from_down;
    }

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_DOWN);
    if (nb >= 0) {
      gf->send_to_down = gfields[nb].recv_from_up;
    }

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_LEFT);
    if (nb >= 0) {
      gf->send_to_left = gfields[nb].recv_from_right;
    }

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_RIGHT);
    if (nb >= 0) {
      gf->send_to_right = gfields[nb].recv_from_left;
    }
  }
}

void group_field_fill(GroupField *gf, double val) {
  if (!gf || !gf->u_prev) return;
  size_t n = (size_t)gf->ny_padded * (size_t)gf->nx_padded;

  if (val == 0.0) {
    memset(gf->u_prev, 0, gf->plane_bytes);
    memset(gf->u_curr, 0, gf->plane_bytes);
    memset(gf->u_next, 0, gf->plane_bytes);
    return;
  }

  for (size_t i = 0; i < n; i++) {
    gf->u_prev[i] = val;
    gf->u_curr[i] = val;
    gf->u_next[i] = val;
  }
}

void group_field_swap(GroupField *gf) {
  if (!gf) return;
  double *tmp = gf->u_prev;
  gf->u_prev = gf->u_curr;
  gf->u_curr = gf->u_next;
  gf->u_next = tmp;
}

void group_field_dump(const GroupField *gf, int gid, FILE *out) {
  if (!gf || !out) return;
  fprintf(out,
    "[GroupField %d] dims=(%d,%d)+2*halo  plane=%.3f MiB  "
    "numa_node=%d  energy=%.6f  "
    "u_prev=%p u_curr=%p u_next=%p\n",
    gid, gf->ny_padded, gf->nx_padded,
    (double)gf->plane_bytes / (1024.0 * 1024.0),
    gf->numa_node, gf->group_energy,
    (void*)gf->u_prev, (void*)gf->u_curr, (void*)gf->u_next);
}
