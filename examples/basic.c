#include "ansa/fxp_fndsa.h"
#include <stdio.h>
#include <stdlib.h>
int main(void) {
    unsigned n = 512;
    size_t pl=fxp_fndsa_public_key_bytes(n), sl=fxp_fndsa_secret_key_bytes(n);
    size_t cap=fxp_fndsa_signature_bytes(n), written=0;
    unsigned char *pk=malloc(pl), *sk=malloc(sl), *sig=malloc(cap);
    fxp_fndsa_prepared *key=NULL;
    const unsigned char message[]="Fixed-point FN-DSA";
    int status=FXP_FNDSA_MEMORY;
    if (!pk || !sk || !sig) goto done;
    status=fxp_fndsa_keygen(n,NULL,pk,pl,sk,sl);
    if (status) goto done;
    status=fxp_fndsa_prepare(&key,sk,sl);
    if (status) goto done;
    status=fxp_fndsa_sign(key,message,sizeof message-1,NULL,sig,cap,&written);
    if (status) goto done;
    status=fxp_fndsa_verify(pk,pl,message,sizeof message-1,sig,written);
    if (!status) printf("Verified %zu-byte signature\n",written);
done:
    fxp_fndsa_free(key);
    if (sk) { volatile unsigned char *p=sk; for(size_t i=0;i<sl;i++) p[i]=0; }
    free(sk);free(pk);free(sig);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
