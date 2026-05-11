#include "mythread/mythread_config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

static int tests_run = 0, tests_passed = 0, tests_failed = 0;

#define TEST(n)  do { tests_run++; printf("  RUN  %s ... ", n); fflush(stdout); } while(0)
#define PASS()   do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(m)  do { tests_failed++; printf("FAIL: %s\n", m); return; } while(0)
#define CHK(c,m) do { if (!(c)) { FAIL(m); } } while(0)

static char *write_temp_cfg(const char *content) {
  char path[] = "/tmp/mythread_test_cfg_XXXXXX";
  int fd = mkstemp(path);
  if (fd < 0) return NULL;
  size_t len = strlen(content);
  if (write(fd, content, len) != (ssize_t)len) { close(fd); return NULL; }
  close(fd);
  return strdup(path);
}

/* TC-C01: load valid config, query in section */
static void test_cfg_load_valid(void) {
  TEST("cfg-load-valid");
  const char *content =
    "[server]\n"
    "host = 127.0.0.1\n"
    "port = 8080\n"
    "\n"
    "[worker]\n"
    "threads = 4\n"
    "timeout = 30.5\n";
  char *path = write_temp_cfg(content);
  mythread_cfg *cfg = mythread_cfg_load(path);
  CHK(cfg != NULL, "cfg_load returned NULL");

  const char *host = mythread_cfg_get(cfg, "server", "host", "default");
  CHK(host != NULL && strcmp(host, "127.0.0.1") == 0, "server.host mismatch");

  int port = mythread_cfg_get_int(cfg, "server", "port", 0);
  CHK(port == 8080, "server.port mismatch");

  int threads = mythread_cfg_get_int(cfg, "worker", "threads", 0);
  CHK(threads == 4, "worker.threads mismatch");

  double timeout = mythread_cfg_get_double(cfg, "worker", "timeout", 0.0);
  CHK(timeout == 30.5, "worker.timeout mismatch");

  mythread_cfg_free(cfg);
  unlink(path);
  free(path);
  PASS();
}

/* TC-C02: missing file returns NULL */
static void test_cfg_load_missing(void) {
  TEST("cfg-load-missing-file");
  mythread_cfg *cfg = mythread_cfg_load("/tmp/__no_such_file__.cfg");
  CHK(cfg == NULL, "missing file should return NULL");
  PASS();
}

/* TC-C03: default fallback for missing key */
static void test_cfg_get_default(void) {
  TEST("cfg-get-default");
  const char *content = "[s]\nk = v\n";
  char *path = write_temp_cfg(content);
  mythread_cfg *cfg = mythread_cfg_load(path);

  CHK(strcmp(mythread_cfg_get(cfg, "s", "missing", "def"), "def") == 0,
      "default string mismatch");
  CHK(mythread_cfg_get_int(cfg, "s", "missing", 42) == 42,
      "default int mismatch");
  CHK(mythread_cfg_get_double(cfg, "s", "missing", 3.14) == 3.14,
      "default double mismatch");

  mythread_cfg_free(cfg);
  unlink(path);
  free(path);
  PASS();
}

/* TC-C04: NULL section treated as global */
static void test_cfg_null_section(void) {
  TEST("cfg-null-section");
  const char *content = "key1 = val1\n[sec]\nkey2 = val2\n";
  char *path = write_temp_cfg(content);
  mythread_cfg *cfg = mythread_cfg_load(path);

  /* key before any [section] is in "" */
  CHK(strcmp(mythread_cfg_get(cfg, NULL, "key1", "def"), "val1") == 0,
      "null section key1 mismatch");
  /* key in section */
  CHK(strcmp(mythread_cfg_get(cfg, "sec", "key2", "def"), "val2") == 0,
      "sec.key2 mismatch");

  mythread_cfg_free(cfg);
  unlink(path);
  free(path);
  PASS();
}

/* TC-C05: NULL cfg returns default */
static void test_cfg_null_cfg(void) {
  TEST("cfg-null-cfg");
  CHK(strcmp(mythread_cfg_get(NULL, "s", "k", "x"), "x") == 0,
      "NULL cfg get string failed");
  CHK(mythread_cfg_get_int(NULL, "s", "k", 7) == 7,
      "NULL cfg get int failed");
  CHK(mythread_cfg_get_double(NULL, "s", "k", 2.5) == 2.5,
      "NULL cfg get double failed");
  mythread_cfg_free(NULL); /* should not crash */
  PASS();
}

/* TC-C06: comments and blank lines ignored */
static void test_cfg_comments(void) {
  TEST("cfg-comments-blanks");
  const char *content =
    "# this is a comment\n"
    "   # indented comment\n"
    "\n"
    "[sec]\n"
    "  val  =   123  \n"
    "  # inline comment on its own line\n";
  char *path = write_temp_cfg(content);
  mythread_cfg *cfg = mythread_cfg_load(path);

  CHK(mythread_cfg_get_int(cfg, "sec", "val", 0) == 123,
      "value after comments mismatch");

  mythread_cfg_free(cfg);
  unlink(path);
  free(path);
  PASS();
}

/* TC-C07: env_get returns env var or default */
static void test_env_get(void) {
  TEST("env-get");
  setenv("MYTHREAD_TEST_ENV", "hello", 1);
  CHK(strcmp(mythread_env_get("MYTHREAD_TEST_ENV", "def"), "hello") == 0,
      "env_get mismatch");
  CHK(strcmp(mythread_env_get("MYTHREAD_NO_SUCH_VAR", "fallback"), "fallback") == 0,
      "env_get default mismatch");
  unsetenv("MYTHREAD_TEST_ENV");
  PASS();
}

/* TC-C08: env_get_int */
static void test_env_get_int(void) {
  TEST("env-get-int");
  setenv("MYTHREAD_TEST_NUM", "99", 1);
  CHK(mythread_env_get_int("MYTHREAD_TEST_NUM", 0) == 99,
      "env_get_int mismatch");
  CHK(mythread_env_get_int("MYTHREAD_NO_SUCH_NUM", -1) == -1,
      "env_get_int default mismatch");
  unsetenv("MYTHREAD_TEST_NUM");
  PASS();
}

/* TC-C09: invalid int falls back to default */
static void test_cfg_invalid_int(void) {
  TEST("cfg-invalid-int");
  const char *content = "[s]\nbad = not_a_number\n";
  char *path = write_temp_cfg(content);
  mythread_cfg *cfg = mythread_cfg_load(path);

  CHK(mythread_cfg_get_int(cfg, "s", "bad", 77) == 77,
      "invalid int should return default");

  mythread_cfg_free(cfg);
  unlink(path);
  free(path);
  PASS();
}

/* TC-C10: bad env int falls back to default */
static void test_env_invalid_int(void) {
  TEST("env-invalid-int");
  setenv("MYTHREAD_BAD_NUM", "not_a_number", 1);
  CHK(mythread_env_get_int("MYTHREAD_BAD_NUM", 55) == 55,
      "invalid env int should return default");
  unsetenv("MYTHREAD_BAD_NUM");
  PASS();
}

/* TC-C11: empty file returns valid cfg with 0 entries */
static void test_cfg_empty_file(void) {
  TEST("cfg-empty-file");
  char *path = write_temp_cfg("");
  mythread_cfg *cfg = mythread_cfg_load(path);
  CHK(cfg != NULL, "empty file should return valid cfg");

  /* queries should all default */
  CHK(strcmp(mythread_cfg_get(cfg, "s", "k", "d"), "d") == 0,
      "empty cfg get failed");

  mythread_cfg_free(cfg);
  unlink(path);
  free(path);
  PASS();
}

int main(void) {
  printf("\n=== mythread config tests ===\n\n");

  test_cfg_load_valid();
  test_cfg_load_missing();
  test_cfg_get_default();
  test_cfg_null_section();
  test_cfg_null_cfg();
  test_cfg_comments();
  test_env_get();
  test_env_get_int();
  test_cfg_invalid_int();
  test_env_invalid_int();
  test_cfg_empty_file();

  printf("\n=== Results: %d run, %d passed, %d failed ===\n",
         tests_run, tests_passed, tests_failed);
  return tests_failed > 0 ? 1 : 0;
}
