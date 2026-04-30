#ifndef MTHREAD_DECOMP_H_INCLUDED
#define MTHREAD_DECOMP_H_INCLUDED
#include "mythread_types.h"

/*
 * 通用域分解接口
 *
 * 将二维矩形域按配置的策略分配给线程组和组内 worker，
 * 自动建立组间邻居拓扑，供 stencil/PDE 类应用直接使用。
 *
 * 使用流程:
 *   mythread_decomp *dc = mythread_decomp_create(nx, ny, halo,
 *       ngroups, nworkers_per_group, DECOMP_Y_ONLY);
 *   // 遍历 dc->group_tiles[gid] 获取每组子域
 *   // 通过 dc->group_neighbors[gid*4 + dir] 查询邻居
 *   // 遍历 dc->worker_tiles[tid] 获取 worker 子域
 *   mythread_decomp_free(dc);
 */

/* ── 分解策略 ── */
typedef enum {
  MYTHREAD_DECOMP_Y_ONLY = 0,   /* 仅 Y 方向切分组，每组覆盖全部 X */
  MYTHREAD_DECOMP_XY_2D  = 1    /* X+Y 都切分，组排成 2D 网格 */
} mythread_decomp_policy;

/* ── 邻居方向 ── */
typedef enum {
  MYTHREAD_NEIGHBOR_UP    = 0,
  MYTHREAD_NEIGHBOR_DOWN  = 1,
  MYTHREAD_NEIGHBOR_LEFT  = 2,
  MYTHREAD_NEIGHBOR_RIGHT = 3,
  MYTHREAD_NEIGHBOR_COUNT = 4
} mythread_neighbor_dir;

/* ── 矩形 tile 描述 ── */
typedef struct {
  int x_begin, x_end;      /* 内部区域（不含 halo），左闭右开 [begin, end) */
  int y_begin, y_end;
  int nx, ny;              /* 内部区域尺寸 */
} mythread_tile;

/* ── 域分解结果 ── */
typedef struct {
  /* 输入参数 */
  int domain_nx, domain_ny;       /* 全域尺寸（内部点） */
  int halo;                       /* halo 宽度 */
  int n_groups;                   /* 组数 */
  int n_workers_per_group;        /* 每组 worker 数（不含 GMT） */
  mythread_decomp_policy policy;  /* 分解策略 */

  /* 组级分解 */
  int gx, gy;                     /* 组的 2D 网格排布（gx * gy = n_groups） */
  mythread_tile *group_tiles;     /* [n_groups]，每组的子域 tile */

  /* 组间邻居拓扑 */
  int *group_neighbors;           /* [n_groups * 4]，flat: [gid*4 + dir] = neighbor_gid 或 -1 */

  /* 组内 worker 分解 */
  mythread_tile *worker_tiles;    /* [n_groups * n_workers_per_group] */
                                  /* flat: worker_tiles[gid * n_workers + tid] */
} mythread_decomp;

/* ── API ── */

/*
 * 创建域分解。
 * nworkers_per_group 可以为 0，此时不计算 worker_tiles。
 * 返回 NULL 表示参数无效。
 */
mythread_decomp *mythread_decomp_create(int domain_nx, int domain_ny, int halo,
                                        int n_groups, int n_workers_per_group,
                                        mythread_decomp_policy policy);

/* 释放域分解 */
void mythread_decomp_free(mythread_decomp *dc);

/*
 * 查询某个组在指定方向上的邻居组 id。
 * 返回 -1 表示该方向无邻居（域边界或 MPI 边界）。
 */
static inline int mythread_decomp_neighbor(const mythread_decomp *dc,
                                           int gid, mythread_neighbor_dir dir) {
  return dc->group_neighbors[gid * MYTHREAD_NEIGHBOR_COUNT + (int)dir];
}

/*
 * 判断某个组在指定方向上是否为域边界。
 * 即：该方向无组内邻居，需要 MPI 通信。
 */
static inline int mythread_decomp_is_domain_boundary(const mythread_decomp *dc,
                                                     int gid,
                                                     mythread_neighbor_dir dir) {
  return mythread_decomp_neighbor(dc, gid, dir) < 0;
}

/*
 * 获取 worker tile。
 * gid: 组 id (0..n_groups-1)
 * wid: 组内 worker 索引 (0..n_workers_per_group-1)
 */
static inline const mythread_tile *mythread_decomp_worker_tile(
    const mythread_decomp *dc, int gid, int wid) {
  return &dc->worker_tiles[gid * dc->n_workers_per_group + wid];
}

/*
 * 打印分解信息（调试用）。
 */
void mythread_decomp_dump(const mythread_decomp *dc, FILE *out);

#endif /* MTHREAD_DECOMP_H_INCLUDED */
