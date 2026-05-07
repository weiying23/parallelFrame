#include "mythread_config.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE   512
#define MAX_ENTRIES 256

typedef struct {
  char section[128];
  char key[128];
  char val[256];
} cfg_entry;

struct mythread_cfg {
  cfg_entry entries[MAX_ENTRIES];
  int n_entries;
};

mythread_cfg *mythread_cfg_load(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) return NULL;

  mythread_cfg *cfg = (mythread_cfg*)calloc(1, sizeof(*cfg));
  if (!cfg) { fclose(f); return NULL; }

  char current_section[128] = "";
  char line[MAX_LINE];

  while (fgets(line, sizeof(line), f)) {
    size_t len = strlen(line);
    while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
      line[--len] = '\0';

    char *s = line;
    while (*s && isspace((unsigned char)*s)) s++;
    if (*s == '\0' || *s == '#') continue;

    /* [section] */
    if (*s == '[') {
      char *end = strchr(s, ']');
      if (end) {
        size_t n = (size_t)(end - s - 1);
        if (n >= sizeof(current_section)) n = sizeof(current_section) - 1;
        memcpy(current_section, s + 1, n);
        current_section[n] = '\0';
      }
      continue;
    }

    /* key = value */
    char *eq = strchr(s, '=');
    if (!eq || cfg->n_entries >= MAX_ENTRIES) continue;

    char *ks = s, *ke = eq;
    while (ke > ks && isspace((unsigned char)ke[-1])) ke--;
    size_t klen = (size_t)(ke - ks);
    if (klen >= 128) klen = 127;

    char *vs = eq + 1;
    while (*vs && isspace((unsigned char)*vs)) vs++;
    char *ve = vs + strlen(vs);
    while (ve > vs && isspace((unsigned char)ve[-1])) ve--;
    size_t vlen = (size_t)(ve - vs);
    if (vlen >= 256) vlen = 255;

    if (klen == 0) continue;

    cfg_entry *e = &cfg->entries[cfg->n_entries++];
    memcpy(e->section, current_section, sizeof(e->section) - 1);
    e->section[sizeof(e->section) - 1] = '\0';
    memcpy(e->key, ks, klen); e->key[klen] = '\0';
    memcpy(e->val, vs, vlen); e->val[vlen] = '\0';
  }

  fclose(f);
  return cfg;
}

void mythread_cfg_free(mythread_cfg *cfg) {
  free(cfg);
}

static const char *cfg_find(const mythread_cfg *cfg,
                            const char *section, const char *key) {
  if (!cfg || !key) return NULL;
  if (!section) section = "";
  for (int i = 0; i < cfg->n_entries; i++) {
    if (strcmp(cfg->entries[i].key, key) == 0 &&
        strcmp(cfg->entries[i].section, section) == 0)
      return cfg->entries[i].val;
  }
  return NULL;
}

const char *mythread_cfg_get(const mythread_cfg *cfg,
                             const char *section, const char *key,
                             const char *defval) {
  const char *v = cfg_find(cfg, section, key);
  return v ? v : defval;
}

int mythread_cfg_get_int(const mythread_cfg *cfg,
                         const char *section, const char *key,
                         int defval) {
  const char *v = cfg_find(cfg, section, key);
  if (!v) return defval;
  char *end;
  long n = strtol(v, &end, 0);
  return (end > v) ? (int)n : defval;
}

double mythread_cfg_get_double(const mythread_cfg *cfg,
                               const char *section, const char *key,
                               double defval) {
  const char *v = cfg_find(cfg, section, key);
  if (!v) return defval;
  char *end;
  double d = strtod(v, &end);
  return (end > v) ? d : defval;
}

const char *mythread_env_get(const char *name, const char *defval) {
  const char *v = getenv(name);
  return v ? v : defval;
}

int mythread_env_get_int(const char *name, int defval) {
  const char *v = getenv(name);
  if (!v) return defval;
  char *end;
  long n = strtol(v, &end, 0);
  return (end > v) ? (int)n : defval;
}
