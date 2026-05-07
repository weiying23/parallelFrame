#include "mythread_halo.h"
#include <mpi.h>
#include <string.h>

/* ═══════════════════════════════════════════════════════════════
 *  组内进程 halo 交换（纯本地 memcpy）
 * ═══════════════════════════════════════════════════════════════ */

void mythread_halo_exchange_intra(GroupField *gfields, const mythread_decomp *dc) {
  if (!gfields || !dc) return;

  const int halo = dc->halo;
  const int ngrp = dc->n_groups;

  /* Phase 1: 发送 — 拷贝源组内边界到目标组 recv 缓冲 */
  for (int g = 0; g < ngrp; g++) {
    GroupField *gf = &gfields[g];
    if (!gf->u_curr) continue;

    const mythread_tile *tile = &dc->group_tiles[g];
    int nx = tile->nx, ny = tile->ny;
    int nb;

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_UP);
    if (nb >= 0 && gfields[nb].recv_from_down) {
      /* g 的末行 → 上方邻居的 DOWN halo（邻居接收来自下方的数据） */
      memcpy(gfields[nb].recv_from_down,
             &gf->u_curr[ny * gf->stride + halo],
             (size_t)nx * sizeof(double));
    }

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_DOWN);
    if (nb >= 0 && gfields[nb].recv_from_up) {
      /* g 的首行 → 下方邻居的 UP halo（邻居接收来自上方的数据） */
      memcpy(gfields[nb].recv_from_up,
             &gf->u_curr[halo * gf->stride + halo],
             (size_t)nx * sizeof(double));
    }

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_LEFT);
    if (nb >= 0 && gfields[nb].recv_from_right) {
      for (int row = 0; row < ny; row++)
        gfields[nb].recv_from_right[row] = gf->u_curr[(halo + row) * gf->stride + halo];
    }

    nb = mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_RIGHT);
    if (nb >= 0 && gfields[nb].recv_from_left) {
      int src_col = halo + nx - 1;
      for (int row = 0; row < ny; row++)
        gfields[nb].recv_from_left[row] = gf->u_curr[(halo + row) * gf->stride + src_col];
    }
  }

  /* Phase 2: 接收 — 从 recv 缓冲应用到 halo */
  for (int g = 0; g < ngrp; g++) {
    GroupField *gf = &gfields[g];
    if (!gf->u_curr) continue;

    const mythread_tile *tile = &dc->group_tiles[g];
    int nx = tile->nx, ny = tile->ny;

    /*
     * GroupField 坐标: y=0 映射到最小全局 Y(下方), y=ny+HALO 映射到最大全局 Y(上方)
     * → recv_from_up(来自上方) 写入 row ny+HALO, recv_from_down(来自下方) 写入 row 0
     */
    if (gf->recv_from_up) {
      memcpy(&gf->u_curr[(ny + halo) * gf->stride + halo], gf->recv_from_up,
             (size_t)nx * sizeof(double));
    }
    if (gf->recv_from_down) {
      memcpy(&gf->u_curr[0 * gf->stride + halo], gf->recv_from_down,
             (size_t)nx * sizeof(double));
    }
    if (gf->recv_from_left) {
      for (int row = 0; row < ny; row++)
        gf->u_curr[(halo + row) * gf->stride + 0] = gf->recv_from_left[row];
    }
    if (gf->recv_from_right) {
      int dst_col = halo + nx;
      for (int row = 0; row < ny; row++)
        gf->u_curr[(halo + row) * gf->stride + dst_col] = gf->recv_from_right[row];
    }
  }
}

/* ═══════════════════════════════════════════════════════════════
 *  域边界 MPI halo 交换
 * ═══════════════════════════════════════════════════════════════ */

void mythread_halo_exchange_mpi(GroupField *gfields, const mythread_decomp *dc,
                                const mythread_mpi_ctx *ctx, uintptr_t comm_arg) {
  if (!gfields || !dc || !ctx) return;

  MPI_Comm comm = (MPI_Comm)comm_arg;
  const int halo = dc->halo;
  const int ngrp = dc->n_groups;
  const int tag_base = ctx->mpi_tag_base > 0 ? ctx->mpi_tag_base : 100;

  /* Y 上边界：没有 UP 邻居的组 */
  if (ctx->mpirank_up >= 0) {
    for (int g = 0; g < ngrp; g++) {
      GroupField *gf = &gfields[g];
      if (!gf->u_curr || !gf->send_up || !gf->recv_up) continue;

      const mythread_tile *tile = &dc->group_tiles[g];
      int nx = tile->nx, ny = tile->ny;

      memcpy(gf->send_up, &gf->u_curr[ny * gf->stride + halo],
             (size_t)nx * sizeof(double));

      MPI_Request req;
      MPI_Isend(gf->send_up, nx, MPI_DOUBLE, ctx->mpirank_up, tag_base,
                comm, &req);
      MPI_Recv(gf->recv_up, nx, MPI_DOUBLE, ctx->mpirank_up, tag_base + 1,
               comm, MPI_STATUS_IGNORE);
      MPI_Wait(&req, MPI_STATUS_IGNORE);

      memcpy(&gf->u_curr[(ny + halo) * gf->stride + halo], gf->recv_up,
             (size_t)nx * sizeof(double));
      break;
    }
  }

  /* Y 下边界：没有 DOWN 邻居的组 */
  if (ctx->mpirank_down >= 0) {
    for (int g = 0; g < ngrp; g++) {
      GroupField *gf = &gfields[g];
      if (!gf->u_curr || !gf->send_down || !gf->recv_down) continue;

      const mythread_tile *tile = &dc->group_tiles[g];
      int nx = tile->nx;

      memcpy(gf->send_down, &gf->u_curr[halo * gf->stride + halo],
             (size_t)nx * sizeof(double));

      MPI_Request req;
      MPI_Isend(gf->send_down, nx, MPI_DOUBLE, ctx->mpirank_down, tag_base + 1,
                comm, &req);
      MPI_Recv(gf->recv_down, nx, MPI_DOUBLE, ctx->mpirank_down, tag_base,
               comm, MPI_STATUS_IGNORE);
      MPI_Wait(&req, MPI_STATUS_IGNORE);

      memcpy(&gf->u_curr[0 * gf->stride + halo], gf->recv_down,
             (size_t)nx * sizeof(double));
      break;
    }
  }

  /* X 左边界：所有没有 LEFT 邻居的组 */
  if (ctx->mpirank_left >= 0) {
    for (int g = 0; g < ngrp; g++) {
      GroupField *gf = &gfields[g];
      if (!gf->u_curr || !gf->send_left || !gf->recv_left) continue;
      if (mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_LEFT) >= 0) continue;

      const mythread_tile *tile = &dc->group_tiles[g];
      int ny = tile->ny;

      for (int row = 0; row < ny; row++)
        gf->send_left[row] = gf->u_curr[(halo + row) * gf->stride + halo];

      MPI_Request req;
      MPI_Isend(gf->send_left, ny, MPI_DOUBLE, ctx->mpirank_left, tag_base + 2,
                comm, &req);
      MPI_Recv(gf->recv_left, ny, MPI_DOUBLE, ctx->mpirank_left, tag_base + 3,
               comm, MPI_STATUS_IGNORE);
      MPI_Wait(&req, MPI_STATUS_IGNORE);

      for (int row = 0; row < ny; row++)
        gf->u_curr[(halo + row) * gf->stride + 0] = gf->recv_left[row];
    }
  }

  /* X 右边界 */
  if (ctx->mpirank_right >= 0) {
    for (int g = 0; g < ngrp; g++) {
      GroupField *gf = &gfields[g];
      if (!gf->u_curr || !gf->send_right || !gf->recv_right) continue;
      if (mythread_decomp_neighbor(dc, g, MYTHREAD_NEIGHBOR_RIGHT) >= 0) continue;

      const mythread_tile *tile = &dc->group_tiles[g];
      int nx = tile->nx, ny = tile->ny;
      int src_col = halo + nx - 1;

      for (int row = 0; row < ny; row++)
        gf->send_right[row] = gf->u_curr[(halo + row) * gf->stride + src_col];

      MPI_Request req;
      MPI_Isend(gf->send_right, ny, MPI_DOUBLE, ctx->mpirank_right, tag_base + 3,
                comm, &req);
      MPI_Recv(gf->recv_right, ny, MPI_DOUBLE, ctx->mpirank_right, tag_base + 2,
               comm, MPI_STATUS_IGNORE);
      MPI_Wait(&req, MPI_STATUS_IGNORE);

      int dst_col = halo + nx;
      for (int row = 0; row < ny; row++)
        gf->u_curr[(halo + row) * gf->stride + dst_col] = gf->recv_right[row];
    }
  }
}
