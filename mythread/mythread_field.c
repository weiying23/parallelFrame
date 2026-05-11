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

/*
 * sched_getcpu() 声明在 <sched.h>，需 _GNU_SOURCE。
 * <numa.h> 可能间接引入，但不应依赖于此。
 */
#if HAS_LIBNUMA
# if __has_include(<sched.h>)
#  include <sched.h>
# endif
#endif

/* ── 平台无关的分配包装 ── */
static int g_field_alloc_warned = 0;

static void *field_alloc(size_t bytes, int node, int *fell_back) {
  if (bytes == 0) return NULL;
  if (fell_back) *fell_back = 0;
#if HAS_LIBNUMA
  if (node >= 0 && numa_available() >= 0) {
    void *p = numa_alloc_onnode(bytes, node);
    if (p) { memset(p, 0, bytes); return p; }
    if (!g_field_alloc_warned) {
      g_field_alloc_warned = 1;
      fprintf(stderr,
        "[mythread] numa_alloc_onnode(%zu bytes, node %d) failed — "
        "falling back to malloc (performance may degrade)\n", bytes, node);
    }
    if (fell_back) *fell_back = 1;
  }
#endif
  void *p = hmalloc(bytes);
  if (p) memset(p, 0, bytes);
  return p;
}

static void field_free(void *p, size_t bytes, int numa_ok) {
  if (!p) return;
#if HAS_LIBNUMA
  if (numa_ok && numa_available() >= 0) {
    numa_free(p, bytes);
    return;
  }
#endif
  (void)bytes;
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
  int halo = dc->halo, nx = tile->nx, ny = tile->ny;
  int node = get_numa_node(), fb = 0, fell_back = 0;

  memset(gf, 0, sizeof(*gf));

  gf->ny_padded   = ny + 2 * halo;
  gf->nx_padded   = nx + 2 * halo;
  gf->stride      = gf->nx_padded;
  gf->plane_bytes = (size_t)gf->ny_padded * (size_t)gf->nx_padded * sizeof(double);
  gf->numa_node   = node;
  gf->numa_ok     = 1;
  gf->group_energy = 0.0;

  gf->u_prev = (double*)field_alloc(gf->plane_bytes, node, &fb); if (fb) fell_back=1;
  gf->u_curr = (double*)field_alloc(gf->plane_bytes, node, &fb); if (fb) fell_back=1;
  gf->u_next = (double*)field_alloc(gf->plane_bytes, node, &fb); if (fb) fell_back=1;
  if (!gf->u_prev || !gf->u_curr || !gf->u_next) goto fail;

#define ALLOC_BUF(ptr, sz) do { \
    (ptr) = (double*)field_alloc((sz), node, &fb); if (fb) fell_back=1; \
    if (!(ptr)) goto fail; \
  } while(0)

  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_UP) >= 0)
    ALLOC_BUF(gf->recv_from_up, (size_t)nx * sizeof(double));
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_DOWN) >= 0)
    ALLOC_BUF(gf->recv_from_down, (size_t)nx * sizeof(double));
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_LEFT) >= 0)
    ALLOC_BUF(gf->recv_from_left, (size_t)ny * sizeof(double));
  if (mythread_decomp_neighbor(dc, gid, MYTHREAD_NEIGHBOR_RIGHT) >= 0)
    ALLOC_BUF(gf->recv_from_right, (size_t)ny * sizeof(double));

  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_UP)) {
    ALLOC_BUF(gf->send_up, (size_t)nx * sizeof(double));
    ALLOC_BUF(gf->recv_up, (size_t)nx * sizeof(double));
  }
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_DOWN)) {
    ALLOC_BUF(gf->send_down, (size_t)nx * sizeof(double));
    ALLOC_BUF(gf->recv_down, (size_t)nx * sizeof(double));
  }
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_LEFT)) {
    ALLOC_BUF(gf->send_left, (size_t)ny * sizeof(double));
    ALLOC_BUF(gf->recv_left, (size_t)ny * sizeof(double));
  }
  if (mythread_decomp_is_domain_boundary(dc, gid, MYTHREAD_NEIGHBOR_RIGHT)) {
    ALLOC_BUF(gf->send_right, (size_t)ny * sizeof(double));
    ALLOC_BUF(gf->recv_right, (size_t)ny * sizeof(double));
  }
#undef ALLOC_BUF

  if (fell_back) gf->numa_ok = 0;
  return 0;

fail:
  group_free_field(gf);
  return -1;
}

void group_free_field(GroupField *gf) {
  if (!gf) return;
  int nok = gf->numa_ok;

  field_free(gf->u_prev,  gf->plane_bytes, nok);
  field_free(gf->u_curr,  gf->plane_bytes, nok);
  field_free(gf->u_next,  gf->plane_bytes, nok);

  field_free(gf->recv_from_up,    0, nok);
  field_free(gf->recv_from_down,  0, nok);
  field_free(gf->recv_from_left,  0, nok);
  field_free(gf->recv_from_right, 0, nok);

  field_free(gf->send_up,    0, nok);
  field_free(gf->recv_up,    0, nok);
  field_free(gf->send_down,  0, nok);
  field_free(gf->recv_down,  0, nok);
  field_free(gf->send_left,  0, nok);
  field_free(gf->recv_left,  0, nok);
  field_free(gf->send_right, 0, nok);
  field_free(gf->recv_right, 0, nok);

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
