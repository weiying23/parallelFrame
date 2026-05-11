#ifndef MTHREAD_FIELD_H_INCLUDED
#define MTHREAD_FIELD_H_INCLUDED
#include "mythread_decomp.h"

/*
 * NUMA 感知的线程组独立内存分配
 *
 * 每个线程组拥有一份独立的波场平面和 halo 缓冲，
 * 物理内存分配在该组绑定的 NUMA 节点上。
 *
 * 使用流程:
 *   // 1. GMT 线程中（bindcpu 已生效）
 *   group_alloc_field(&g_gfields[gid], gid, g_decomp);
 *
 *   // 2. 访问元素
 *   double v = gf->u_curr[GFIDX(gf, y, x)];
 *
 *   // 3. 时间步推进
 *   group_field_swap(gf);
 *
 *   // 4. 释放
 *   group_free_field(gf);
 */

typedef struct GroupField {
  /* ── 波场平面（每组独立分配）── */
  double *u_prev;
  double *u_curr;
  double *u_next;

  /* ── 派生尺寸 ── */
  int ny_padded;               /* tile->ny + 2*halo */
  int nx_padded;               /* tile->nx + 2*halo */
  int stride;                  /* = nx_padded */

  /* ── 完整平面尺寸（字节）── */
  size_t plane_bytes;

  /* ── 组间 halo 缓冲 ── */
  double *send_to_up;          /* 发送给上组：nx 个 double */
  double *recv_from_up;
  double *send_to_down;
  double *recv_from_down;
  double *send_to_left;        /* 发送给左组：ny 个 double（XY_2D） */
  double *recv_from_left;
  double *send_to_right;
  double *recv_from_right;

  /* ── MPI halo 缓冲（仅域边界组分配）── */
  double *send_up,    *recv_up;
  double *send_down,  *recv_down;
  double *send_left,  *recv_left;
  double *send_right, *recv_right;

  /* ── NUMA 信息 ── */
  int numa_node;
  int numa_ok;                 /* 1=NUMA 分配成功, 0=已回退 malloc */

  /* ── 能量 ── */
  double group_energy;
} GroupField;

/*
 * 为组 gid 分配所有波场平面和 halo 缓冲。
 * 必须在 GMT 线程中调用（此时 bindcpu 已生效，NUMA 节点正确）。
 * dc 是已创建的域分解结果。
 * 返回 0 成功，-1 分配失败。
 */
int group_alloc_field(GroupField *gf, int gid, const mythread_decomp *dc);

/*
 * 释放 group_alloc_field 分配的全部内存。
 * 安全处理 NULL 指针和已释放的字段。
 */
void group_free_field(GroupField *gf);

/*
 * 在多个组的 GroupField 间建立 send/recv 缓冲的指针互连。
 * 必须在所有组的 group_alloc_field 完成后、使用 halo 交换前调用。
 */
void group_field_link_buffers(GroupField *gfields, const mythread_decomp *dc);

/*
 * 线性和索引（含 halo，以 HALO 为原点）。
 * 等价于 (y) * (gf)->stride + (x)
 */
#define GFIDX(gf, y, x) ((size_t)(y) * (size_t)(gf)->stride + (size_t)(x))

/*
 * 将三个平面填充为 val（通常用于初始化或清零 halo 区域）。
 */
void group_field_fill(GroupField *gf, double val);

/*
 * 交换指针：prev ← curr, curr ← next, next ← prev（旧）。
 * 用于时间步推进，不移动数据。
 */
void group_field_swap(GroupField *gf);

/*
 * 输出某组的内存信息（调试用）。
 */
void group_field_dump(const GroupField *gf, int gid, FILE *out);

#endif /* MTHREAD_FIELD_H_INCLUDED */
