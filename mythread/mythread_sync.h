#ifndef MTHREAD_SYNC_H_INCLUDED
#define MTHREAD_SYNC_H_INCLUDED
#include "mythread_types.h"

__BEGIN_DECLS
int GetVInt(int volatile *volatile p);

/* SM — Sub ↔ Main */
void sSetState(int state);
void sWaitState(int state);
void sWaitStater(int state);
void mSetSubs(int state);
void mWaitSubs(int state);
void mWaitSubsr(int state);

/* SG — Sub ↔ Group */
void sWaitGrp(int state);
void sWaitGrpr(int state);
void sSetGrp(int state);
void gWaitSubs(int state);
void gWaitSubsr(int state);
void gSetSubs(int state);

/* GM — Group ↔ Main */
void gWaitMain(int state);
void gWaitMainr(int state);
void gSetMain(int state);
void mWaitGrps(int state);
void mWaitGrpsr(int state);
void mSetGrps(int state);
__END_DECLS

/* SM */
#define SSM  sSetState   (RFB)
#define SWM  sWaitState  (RFB)
#define SWMR sWaitStater (RFB)
/* MS */
#define MSS  mSetSubs    (RFB)
#define MWS  mWaitSubs   (RFB)
#define MWSR mWaitSubsr  (RFB)
/* SG */
#define SSG  sSetGrp     (RFB)
#define SWG  sWaitGrp    (RFB)
#define SWGR sWaitGrpr   (RFB)
/* GS */
#define GSS  gSetSubs    (RFB)
#define GWS  gWaitSubs   (RFB)
#define GWSR gWaitSubsr  (RFB)
/* MG */
#define MSG  mSetGrps    (RFB)
#define MWG  mWaitGrps   (RFB)
#define MWGR mWaitGrpsr  (RFB)
/* GM */
#define GSM  gSetMain    (RFB)
#define GWM  gWaitMain   (RFB)
#define GWMR gWaitMainr  (RFB)

#endif
