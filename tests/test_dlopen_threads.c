#include "ansa/fxp_fndsa.h"
#include <dlfcn.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "FAIL %d: %s\n", __LINE__, #x);                                        \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t condition = PTHREAD_COND_INITIALIZER;
static int ready, resume_thread;
static fxp_fndsa_prepared *key;
static __typeof__(fxp_fndsa_sign) *sign_message;
static void *worker(void *unused)
{
    (void)unused;
    uint8_t seed[64] = {1}, signature[666];
    size_t length;
    CHECK(!sign_message(key, NULL, 0, seed, signature, sizeof signature, &length));
    CHECK(length == sizeof signature);
    pthread_mutex_lock(&lock);
    ready = 1;
    pthread_cond_signal(&condition);
    while (!resume_thread)
        pthread_cond_wait(&condition, &lock);
    pthread_mutex_unlock(&lock);
    return NULL; /* TLS destructors must still point into a loaded image. */
}
int main(int argc, char **argv)
{
    CHECK(argc == 2);
    void *library = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
    CHECK(library);
    __typeof__(fxp_fndsa_keygen) *keygen = dlsym(library, "fxp_fndsa_keygen");
    __typeof__(fxp_fndsa_prepare) *prepare = dlsym(library, "fxp_fndsa_prepare");
    __typeof__(fxp_fndsa_free) *free_key = dlsym(library, "fxp_fndsa_free");
    sign_message = dlsym(library, "fxp_fndsa_sign");
    CHECK(keygen && prepare && free_key && sign_message);
    uint8_t seed[64] = {0}, pk[897], sk[1345];
    CHECK(!keygen(512, seed, pk, sizeof pk, sk, sizeof sk));
    CHECK(!prepare(&key, sk, sizeof sk));
    pthread_t thread;
    CHECK(!pthread_create(&thread, NULL, worker, NULL));
    pthread_mutex_lock(&lock);
    while (!ready)
        pthread_cond_wait(&condition, &lock);
    free_key(key);
    CHECK(!dlclose(library));
    resume_thread = 1;
    pthread_cond_signal(&condition);
    pthread_mutex_unlock(&lock);
    CHECK(!pthread_join(thread, NULL));
    puts("PASS: close library handle before worker exits; TLS cleanup remains valid");
    return 0;
}
