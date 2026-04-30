#ifndef MTHREAD_HALO_H_INCLUDED
#define MTHREAD_HALO_H_INCLUDED
#include <stdint.h>
#include "mythread_field.h"
#include "mythread_decomp.h"

/*
 * 组间 + MPI halo 交换
 *
 * Layer 3: 基于 GroupField 和 mythread_decomp 的邻居拓扑，
 * 自动完成组间本地拷贝和域边界 MPI 通信。
 *
 * MPI_Comm 通过 uintptr_t 传递以兼容 MPICH(int) 和 OpenMPI(pointer):
 *   调用方: mythread_halo_exchange_mpi(gfs, dc, &ctx, (uintptr_t)MPI_COMM_WORLD);
 *   实现方: MPI_Comm comm = (MPI_Comm)comm_arg;
 *
 * 使用流程:
 *   group_field_link_buffers(gfields, dc);
 *   mythread_halo_exchange_intra(gfields, dc);
 *   mythread_halo_exchange_mpi(gfields, dc, &mpi_ctx, (uintptr_t)MPI_COMM_WORLD);
 */

typedef struct {
  int mpirank_up;
  int mpirank_down;
  int mpirank_left;
  int mpirank_right;
  int mpi_tag_base;      /* 默认 100 */
} mythread_mpi_ctx;

/* 纯本地 memcpy，不涉及 MPI */
void mythread_halo_exchange_intra(GroupField *gfields, const mythread_decomp *dc);

/* 域边界 MPI 收发，comm_arg 是 (uintptr_t)MPI_Comm */
void mythread_halo_exchange_mpi(GroupField *gfields, const mythread_decomp *dc,
                                const mythread_mpi_ctx *ctx, uintptr_t comm_arg);

/* 完整交换：先 intra 后 mpi */
static inline void mythread_halo_exchange_all(GroupField *gfields,
                                              const mythread_decomp *dc,
                                              const mythread_mpi_ctx *ctx,
                                              uintptr_t comm_arg) {
  mythread_halo_exchange_intra(gfields, dc);
  mythread_halo_exchange_mpi(gfields, dc, ctx, comm_arg);
}

#endif /* MTHREAD_HALO_H_INCLUDED */
