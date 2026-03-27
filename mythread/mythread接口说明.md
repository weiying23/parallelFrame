
# mythread 接口函数文档

## 1. 线程初始化和控制函数

| 函数名 | 参数 | 返回值 | 说明 |
|--------|------|--------|------|
| `InitThreads` | `int mpi_id_, NCorePClu_, NCluPNode_, NCorePGrp_, NThPGrp_, NGrpPProc_, NProcPNode_, int *ManageCoreId_` | `int` | 初始化线程系统，配置线程参数 |
| `StartThreads` | `TFunc tfun` | `void` | 启动所有工作线程，执行指定函数 |
| `EndThreads` | `void` | `void` | 结束所有工作线程 |
| `setthread` | `int NCorePClu_, NThPClu_, NGrpPProc_, NProcPNode_, ManageCoreId_` | `void` | 设置线程配置参数 |
| `initthreads_` | `int *mpi_id_, *NCorePClu_, *NCluPNode_, *NCorePProc_, *NThPGrp_, *NGrpPProc_, *NProcPNode_, *ManageCoreId_, *err` | `void` | Fortran接口：初始化线程 |
| `startthreads_` | `void` | `void` | Fortran接口：启动线程 |
| `endthreads_` | `void` | `void` | Fortran接口：结束线程 |
| `bindthread` | `void` | `void` | 绑定当前线程到指定CPU核心 |
| `bindthread_` | `void` | `void` | Fortran接口：绑定线程 |
| `bindcpu` | `int id` | `int` | 绑定当前线程到指定CPU ID |

---

## 2. 线程信息获取函数

| 函数名 | 参数 | 返回值 | 说明 |
|--------|------|--------|------|
| `gettid` | `void` | `int` | 获取当前线程索引 |
| `gettid_` | `void` | `int` | Fortran接口：获取线程索引 |
| `getnt` | `void` | `int` | 获取当前组线程总数 |
| `getnt_` | `void` | `int` | Fortran接口：获取线程总数 |
| `_threadmain_` | `HTHREADINFO ti` | `void` | 线程主函数入口 |

---

## 3. 同步函数

### 3.1 SubThread <-> MainThread 同步 (SM / MS)

| 函数名 | 参数 | 说明 |
|--------|------|------|
| `sSetState` | `int state` | 设置子线程状态（通知主线程） |
| `sWaitState` | `int state` | 子线程等待主线程状态 |
| `sWaitStater` | `int state` | 子线程等待主线程状态（反向条件） |
| `mSetSubs` | `int state` | 主线程设置子线程状态 |
| `mWaitSubs` | `int state` | 主线程等待所有子线程状态 |
| `mWaitSubsr` | `int state` | 主线程等待所有子线程状态（反向条件） |

**宏定义：**
```c
#define SSM   sSetState(RFB)    // 设置子线程就绪
#define SWM   sWaitState(RFB)   // 等待主线程就绪
#define SWMR  sWaitStater(RFB)  // 等待主线程就绪（反向）
#define MSS   mSetSubs(RFB)     // 主线程设置状态
#define MWS   mWaitSubs(RFB)    // 主线程等待子线程
#define MWSR  mWaitSubsr(RFB)   // 主线程等待子线程（反向）
```

### 3.2 SubThread <-> GroupMain 同步 (SG / GS)

| 函数名 | 参数 | 说明 |
|--------|------|------|
| `sSetGrp` | `int state` | 设置组内子线程状态 |
| `sWaitGrp` | `int state` | 等待组主线程状态 |
| `sWaitGrpr` | `int state` | 等待组主线程状态（反向条件） |
| `gSetSubs` | `int state` | 组主线程设置子线程状态 |
| `gWaitSubs` | `int state` | 组主线程等待所有子线程 |
| `gWaitSubsr` | `int state` | 组主线程等待所有子线程（反向条件） |

**宏定义：**
```c
#define SSG   sSetGrp(RFB)      // 设置组子线程就绪
#define SWG   sWaitGrp(RFB)     // 等待组主线程就绪
#define SWGR  sWaitGrpr(RFB)    // 等待组主线程就绪（反向）
#define GSS   gSetSubs(RFB)     // 组主线程设置状态
#define GWS   gWaitSubs(RFB)    // 组主线程等待子线程
#define GWSR  gWaitSubsr(RFB)   // 组主线程等待子线程（反向）
```

### 3.3 GroupMain <-> MainThread 同步 (GM / MG)

| 函数名 | 参数 | 说明 |
|--------|------|------|
| `gSetMain` | `int state` | 组主线程设置主线程状态 |
| `gWaitMain` | `int state` | 组主线程等待主线程状态 |
| `gWaitMainr` | `int state` | 组主线程等待主线程状态（反向条件） |
| `mSetGrps` | `int state` | 主线程设置所有组状态 |
| `mWaitGrps` | `int state` | 主线程等待所有组主线程 |
| `mWaitGrpsr` | `int state` | 主线程等待所有组主线程（反向条件） |

**宏定义：**
```c
#define GSM   gSetMain(RFB)     // 组主线程设置主线程就绪
#define GWM   gWaitMain(RFB)    // 组主线程等待主线程
#define GWMR  gWaitMainr(RFB)   // 组主线程等待主线程（反向）
#define MSG   mSetGrps(RFB)     // 主线程设置组状态
#define MWG   mWaitGrps(RFB)    // 主线程等待所有组
#define MWGR  mWaitGrpsr(RFB)   // 主线程等待所有组（反向）
```

### 3.4 通用同步

| 函数名 | 参数 | 说明 |
|--------|------|------|
| `MultiThreadSync` | `void` | 多线程同步屏障 |
| `multithreadsync_` | `void` | Fortran接口：多线程同步 |

---

## 4. 性能计时函数 (TSC)

| 函数名 | 参数 | 说明 |
|--------|------|------|
| `tscinit` | `void` | 初始化时间戳计数器 |
| `tscb` | `int id` | 开始计时（标记起点） |
| `tsce` | `int id` | 结束计时（累加时间差） |
| `tsceb` | `int id` | 结束上一计时并开始新计时 |
| `prtsc` | `const char* tag` | 打印时间统计结果 |
| `tscb_` | `int *id_` | Fortran接口：开始计时 |
| `tsce_` | `int *id_` | Fortran接口：结束计时 |
| `tsceb_` | `int *id_` | Fortran接口：结束并开始新计时 |
| `prtsc_` | `void` | Fortran接口：打印时间统计 |

---

## 5. Fortran 接口函数汇总

| C函数 | Fortran接口 |
|-------|-------------|
| `sWaitState` | `swaitstate_` |
| `sWaitStater` | `swaitstater_` |
| `sSetState` | `ssetstate_` |
| `mWaitSubs` | `mwaitsubs_` |
| `mSetSubs` | `msetsubs_` |
| `gWaitSubs` | `gwaitsubs_` |
| `gWaitSubsr` | `gwaitsubsr_` |
| `gSetSubs` | `gsetsubs_` |
| `gWaitMain` | `gwaitmain_` |
| `gWaitMainr` | `gwaitmainr_` |
| `gSetMain` | `gsetmain_` |
| `mWaitGrps` | `mwaitgrps_` |
| `mWaitGrpsr` | `mwaitgrpsr_` |
| `mSetGrps` | `msetgrps_` |
| `sWaitGrp` | `swaitgrp_` |
| `sWaitGrpr` | `swaitgrpr_` |
| `sSetGrp` | `ssetgrp_` |

---

## 6. 其他辅助函数

| 函数名 | 参数 | 返回值 | 说明 |
|--------|------|--------|------|
| `SetLocV` | `int typ, int ind, void* p` | `void` | 设置本地变量指针 |
| `GetLocV` | `int typ, int ind, void* p` | `void*` | 获取本地变量指针 |
| `ntdelay` | `int n` | `void` | 纳秒级延迟 |
| `ntdelay_` | `int n` | `void` | Fortran接口：延迟 |
| `opentf` | `void` | `void` | 打开调试日志文件 |

---

## 7. 全局变量

| 变量名 | 类型 | 说明 |
|--------|------|------|
| `md` | `threadProc` | 全局线程控制结构 |
| `ti` | `THREADINFO*` | 当前线程信息（TLS） |
| `gi` | `threadGroup*` | 当前线程组信息（TLS） |
| `mpi_id` | `int` | MPI进程ID |
| `NCorePClu` | `int` | 每簇CPU核心数 |
| `NCluPNode` | `int` | 每节点簇数 |
| `NCorePGrp` | `int` | 每组CPU核心数 |
| `NGrpPProc` | `int` | 每进程组数 |
| `NThPGrp` | `int` | 每组线程数 |
| `NProcPNode` | `int` | 每节点进程数 |
| `ManageCoreId` | `int` | 管理核心ID |
| `NThreads` | `int` | 总线程数 |
| `ThreadG` | `int` | 分组标志 |

---

## 8. 宏常量

| 宏名 | 默认值 | 说明 |
|------|--------|------|
| `MCOREPC` | 64 | 最大支持的核心数 |
| `MCLUST` | 20 | 最大支持的簇数 |
| `MSB` | 16 | 状态数组大小 |
| `MSBG` | 16 | 组状态数组大小 |

---

## 使用示例

```c
// 初始化线程系统
InitThreads(0, 38, 16, 38, 8, 4, 16, &manageId);

// 启动线程执行计算函数
StartThreads(my_compute_function);

// 在线程函数中使用同步
void my_compute_function() {
    // 子线程通知主线程就绪
    SSM;
    
    // 等待主线程信号
    SWM;
    
    // 执行计算...
}

// 结束线程
EndThreads();
```