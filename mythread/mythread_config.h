#ifndef MTHREAD_CONFIG_H_INCLUDED
#define MTHREAD_CONFIG_H_INCLUDED
#include <stdio.h>

/*
 * 轻量级配置文件解析器。
 *
 * 格式:
 *   # 注释行
 *   key = value
 *   [section]
 *
 * 使用:
 *   mythread_cfg *cfg = mythread_cfg_load("config.cfg");
 *   const char *v = mythread_cfg_get(cfg, "section", "key", "default");
 *   int n = mythread_cfg_get_int(cfg, "section", "key", 100);
 *   mythread_cfg_free(cfg);
 */

typedef struct mythread_cfg mythread_cfg;

/* 加载配置文件，失败返回 NULL */
mythread_cfg *mythread_cfg_load(const char *path);

/* 释放 */
void mythread_cfg_free(mythread_cfg *cfg);

/* 查询字符串值，未找到返回 defval */
const char *mythread_cfg_get(const mythread_cfg *cfg,
                             const char *section, const char *key,
                             const char *defval);

/* 查询整数值 */
int mythread_cfg_get_int(const mythread_cfg *cfg,
                         const char *section, const char *key,
                         int defval);

/* 查询浮点值 */
double mythread_cfg_get_double(const mythread_cfg *cfg,
                               const char *section, const char *key,
                               double defval);

/* 读取环境变量，未设置则返回 defval */
const char *mythread_env_get(const char *name, const char *defval);
int mythread_env_get_int(const char *name, int defval);

#endif /* MTHREAD_CONFIG_H_INCLUDED */
