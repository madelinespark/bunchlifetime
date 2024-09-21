//These routines are executed by the aSub record

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dbDefs.h>
#include <epicsTypes.h>
#include <epicsTime.h>
#include <registryFunction.h>
#include <aSubRecord.h>
#include <epicsExport.h>

#define DEBUG
#ifdef DEBUG
#define DEBUG_PRINT printf
#else
#define DEBUG_PRINT(...)
#endif

#define HARMONIC_NR 328
#define NMAX 1000
#define MAMIN .1

// #define NMAX_LIMIT 2000
//epicsTimeStamp tPrevTime;

// global variables
epicsFloat64 BCM[HARMONIC_NR][NMAX] = {0}; // init everything to 0?
epicsFloat64 Ts[NMAX] = {0};
int iFit0[HARMONIC_NR] = {0};
int iFit1 = 0;
epicsFloat64 Tau[HARMONIC_NR] = {0};
epicsFloat64 I0[HARMONIC_NR] = {0};
unsigned int uiX=0; //index

epicsInt32 TimInjReq0[7];
epicsFloat64 Ts_TimInjReq0=0;
epicsFloat64 dTs_BCM=0;
epicsFloat64 *pdBunchInpWF=NULL;

// helper functions!

void get_row_at(epicsFloat64 buf[], int start, int end, int r, int c, epicsFloat64 arr[][c]) {
    int j = 0;
    for (int i = start; i <= end; i++) { // inclusive of start and end
        buf[j] = arr[r][i];
        j++;
    }
}

void get_column(epicsFloat64 buf[], int nr, int nc, int c, epicsFloat64 arr[][nc]) {
    for (int i = 0; i < nr; i++) {
        buf[i] = arr[i][c];
    }
}

void insert_column(int nr, int nc, int c, epicsFloat64 lst[nr], epicsFloat64 arr[][nc]) {
    for (int i = 0; i < nr; i++) {
        arr[i][c] = lst[i];
    }
}

void low_bunch_current() {
// return index of first arg if smaller than second arg
    epicsFloat64 BCM_col[HARMONIC_NR] = {0};
    
    get_column(BCM_col, HARMONIC_NR, NMAX, uiX, BCM);

    for (int i = 0; i < HARMONIC_NR; i++) {
        if ((BCM_col[i] < MAMIN)) {
            iFit0[i] = uiX + 1;
        }
    }
}

static long bunchLifeTimeProcInit(aSubRecord *pRec)
{
    epicsInt32 *ptr=NULL;
    ptr =  (epicsInt32 *) pRec->c;
    memcpy(TimInjReq0, ptr,sizeof(epicsInt32)*7);
    Ts_TimInjReq0 = *(epicsFloat64 *) pRec->f;
    dTs_BCM = *(epicsFloat64 *) pRec->e;
    pdBunchInpWF = (epicsFloat64 *) pRec->a;

    DEBUG_PRINT("Record %s called bunchLifeTimeProcInit(%p)\n",
                  pRec->name, (void*) pRec);

    /* initialization? */
    insert_column(HARMONIC_NR, NMAX, uiX, pdBunchInpWF, BCM);
    Ts[uiX] = dTs_BCM;


    /* check for low bunch current */
    low_bunch_current();
    return 0;
}

static long bunchLifeTimeProcProcess(aSubRecord *pRec)
{
    epicsInt32 *TimInjReq=NULL;
    epicsFloat64 *pdTimeOutWF=NULL;
    epicsFloat64 *pdAvg = NULL;
//   epicsTimeStamp tCurrentTime;
    epicsFloat64 dTimeDiff=0;
    epicsFloat64 dTs_TimInjReq=0;
    epicsUInt16 sDebug=0;
    epicsInt32 iSrBucketNr=0;
    epicsInt32 iBunchesNr=0;
    epicsInt32 iSeqNr=0;
    epicsFloat64 dNmin=0;
    epicsFloat64 dNmax=0;
    epicsFloat64 dAmin=0;
    int LSFitFlag = 0;
    int iFlag = 0;
    int nFit[HARMONIC_NR] = {0};

    pdBunchInpWF = (epicsFloat64 *) pRec->a;
    sDebug = *(epicsUInt16 *) pRec->b;
    TimInjReq = (epicsInt32 *) pRec->c;
    dTimeDiff = *(epicsFloat64 *) pRec->d;
    dTs_BCM = *(epicsFloat64 *) pRec->e;
    dTs_TimInjReq = *(epicsFloat64 *) pRec->f;
    dNmin = *(epicsFloat64 *) pRec->g;
    dNmax = *(epicsFloat64 *) pRec->h;
    dAmin = *(epicsFloat64 *) pRec->i;

//   tCurrentTime         =  (epicsTimeStamp ) pRec->time;
    pdTimeOutWF = (epicsFloat64 *) pRec->valb;
    pdAvg = (epicsFloat64 *) pRec->valc;
    iSrBucketNr=TimInjReq[0];
    iBunchesNr=TimInjReq[1];
    iSeqNr=TimInjReq[6];


    DEBUG_PRINT("Debug %d\n",sDebug);

    if(sDebug){
        printf("iSrBucketNr=%d\n",iSrBucketNr);
        printf("iBunchesNr=%d\n",iBunchesNr);
        printf("iSeqNr=%d\n",iSeqNr);
        //printf("Time diff [ns]= %lld\n",epicsTimeDiffInNS(&tCurrentTime,&tPrevTime));
        printf("Time diff [s]= %f\n", dTimeDiff);
        printf("Nmin=%f Nmax=%f nAmin=%f\n",dNmin,dNmax,dAmin);
        printf("uiX=%d\n", uiX);
    }
//  tPrevTime.secPastEpoch=tCurrentTime.secPastEpoch;
//   tPrevTime.nsec=tCurrentTime.nsec;
//------  Begining of the Bunch life time calculation algorithm ---

    uiX++;
    /* shift the data if greater than max arr size */
    if (uiX > (dNmax - 1)) {
        uiX = dNmax - 1;

        for (int i = 0; i < HARMONIC_NR; i++) { // shift BCM
            for (int j = 0; j < dNmax - 1; j++) {
                BCM[i][j] = BCM[i][j+1];
            }

            iFit0[i] = ((iFit0[i] - 1) < 0) ? 0 : (iFit0[i] - 1);
        }

        for (int i = 0; i < dNmax - 1; i++) { // shift Ts
            Ts[i] = Ts[i + 1];
        }
        iFit1 = iFit1 - 1;
    }

    insert_column(HARMONIC_NR, dNmax, uiX, pdBunchInpWF, BCM);
    Ts[uiX] = dTs_BCM;
    iFit1 = iFit1 + 1;

    /* check for low bunch current */
    low_bunch_current();

    /* add checks for beam loss */

    /* checks for injection */
    if (TimInjReq[6] != TimInjReq0[6] && TimInjReq[2] == 40 && TimInjReq0[3] == 0) {
        iSrBucketNr=TimInjReq[0]; // target bucket
        iBunchesNr=TimInjReq[1]; // number of bunches (always 1 for now)
        int InjectedBuckets = 0;   // max 4 bunches (?)
        for (int i = 0; i < iBunchesNr; i++) {
            InjectedBuckets = ((4 * i) + iSrBucketNr) - 1;
            iFit0[InjectedBuckets] = uiX + 1; // move iFit0 indexes
        }
    }
        memcpy(TimInjReq0,TimInjReq,sizeof(epicsInt32)*7);
        Ts_TimInjReq0 = dTs_TimInjReq;

    /* LS Fit */
    for (int iBunch = 0; iBunch < HARMONIC_NR; iBunch++) {
        // only fit if there is enough data
        if ((iFit1 - iFit0[iBunch]) > dNmin) {
            int ind = iFit1-iFit0[iBunch] + 1;
            epicsFloat64 d[NMAX] = {0};
            get_row_at(d, iFit0[iBunch], iFit1, iBunch, dNmax, BCM);

            epicsFloat64 t[NMAX] = {0};
            int j = 0;
            for (int i = iFit0[iBunch]; i <= iFit1; i++) {
                t[j] = Ts[i];
                j++;
            }
            
            epicsFloat64 y[NMAX] = {0};
            for (int i = 0; i < ind; i++) {
                y[i] = log(d[i]);
            }
            epicsFloat64 m = ((y[ind - 1] - y[0]) / (t[ind - 1] - t[0]));
            Tau[iBunch] = (-1/m/60/60);

            if (isnan(Tau[iBunch])) {
                printf("Bad lifetime fit for bunch %d, m=%f\n", iBunch, m); // bunch at index iBunch (0-327)
            }
        }
        else {
            I0[iBunch] = NAN;
            Tau[iBunch] = NAN;
        }
    }

    int ii[HARMONIC_NR] = {0};
    int j = 0;
    for (int i = 0; i < HARMONIC_NR; i++) {
        if (!isnan(Tau[i])) {
            if (i <= 190) {
                ii[j] = i; // remove the cam bunch
                j++;
            }
        }
    }

    epicsFloat64 mA[HARMONIC_NR] = {0};
    get_column(mA, HARMONIC_NR, dNmax, (iFit1), BCM);
    
    epicsFloat64 mA_total = 0;
    //printf("%d\n", j);
    epicsFloat64 sum = 0;
    epicsFloat64 lifetime_Avg = 0;
    for (int i = 0; i < j + 1; i++) {
        mA_total += mA[ii[i]];
        sum += (Tau[ii[i]] * mA[ii[i]]);
    }
    lifetime_Avg = sum / mA_total;
    if (isnan(lifetime_Avg)) {
        // setpvonline('SR:BCM:BunchLifetime:Avg', -1)
        *pdAvg = -1;
    }
    else {
        // setpvonline('SR:BCM:BunchLifetime:Avg', lifetime_Avg)
        *pdAvg = lifetime_Avg;
    }

    // look for errors
    if (lifetime_Avg > 20) {
        iFlag = 1;
    }
    else {
        for (int i = 0; i < HARMONIC_NR; i++) {
            if (Tau[i] > 20) {
                iFlag = 1;
                break;
            }
        }
    }

    for (int i = 0; i < HARMONIC_NR; i++) {
        if (isnan(Tau[i])) {
            Tau[i] = -1;
        }
    }

    // setpvonline('SR:BCM:BunchLifetime', Tau')
    for (int i = 0; i < HARMONIC_NR; i++) {
        pdTimeOutWF[i] = Tau[i];
    }


//-----  End of the Bunch life time calculation algorithm ---
    if(sDebug){
        printf("OUT WF [1,2,3]=%f,%f,%f\n",pdTimeOutWF[0],pdTimeOutWF[1],pdTimeOutWF[2]);
        printf("Avg: %f\n\n", *pdAvg);
    }
    
    return 0;
}

// code runs here
epicsRegisterFunction(bunchLifeTimeProcInit);
epicsRegisterFunction(bunchLifeTimeProcProcess);
