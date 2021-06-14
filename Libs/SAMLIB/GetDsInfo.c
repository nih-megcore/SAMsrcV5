// GetDsInfo() -- read ugly clumsy CTF-style data file into
//  beautiful efficient biomag memory structures
//
//  Note: If 'BadOut' is NULL, then do not return 'Bad' array
//
//  Author: SE Robinson (modified for 64-bit by T Holroyd)
//          MEG Core Facility
//          NIMH
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <byteswap.h>
#include <string.h>
#include <errno.h>
#include <geoms.h>
#include <samlib.h>
#include <DataFiles.h>
#include <MEGDefs.h>

#define MAX_MEG4    (2147483640)        // 2^31-8

void GetDsInfo(
    char            *DsName,            // input pathname
    HeaderInfo      *Header,            // header information
    ChannelInfo     **ChannelOut,       // *ChannelOut[M] -- pointer to array of channel information
    EpochInfo       **EpochOut,         // *EpochOut[E] -- pointer to array of epoch information
    unsigned char   ***BadOut,          // *BadOut[E][T] -- pointer to array of sample flags
    int             Xformed             // sensor transform flag: TRUE=head-coil frame, FALSE=sensor frame
) {
    ChannelInfo     *Channel;           // Channel[M] -- channel information
    EpochInfo       *Epoch;             // Epoch[E] -- epoch information
    double          Vc[3];              // Cartesian position variable
    double          *Gain;              // Gain[m] -- array for CTF 'properGain'
    double          ioGain;
    double          properGain;
    double          qGain;
    double          Polarity;           // P-vector
    double          BadStart;           // bad segment starting time
    double          BadEnd;             // bad segment ending time
    double          Dtmp;               // temporary double
    int             c;                  // coil index
    int             i;                  // sample index
    int             j;                  // 'nuther index
    int             e;                  // epoch index
    int             l;                  // string length
    int             m;                  // channel index
    int             n;                  // alternate channel index
    int             t;                  // time index
    int             order;              // balance order index
    int             v;                  // vector index
    int             E;                  // epoch count
    int             M;                  // number channels
    int             T;                  // sample count
    int             NumSens;            // number of primary channels
    int             NumRefs;            // number of reference channels
    int             NumEEGs;            // number of eeg electrodes
    int             NumBad;             // number of bad epochs or channels
    int             NumClasses;         // number of classes
    int             BadEpoch;           // bad epoch number
    int             bs;                 // bad segment starting sample
    int             be;                 // bad segment ending sample
    int             flag;               // label search flag
    int             Itmp;               // temporary integer variable
    size_t          NumBytes;           // number of bytes to read
    unsigned char   **Bad;              // Bad[E][T] -- array of sample flags
    char            StringField[128];   // character array large enought for line
    char            HeaderLabel[8];     // file header label
    char            ChannelLabel[32];   // channel name
    char            DsPath[256];        // dir ds is in
    char            PathName[256];      // full path name
    FILE            *fp;                // generic file pointer
    FILE            *fin;               // input file pointer
    FILE            *feeg;              // optional 'eegl file pointer
    meg41GeneralResRec GeneralResource; // 'meg41' resource header info
    int32_t         Padding;            // 4 byte dummy variable
    char            *run_descriptor;    // should be in header?
    short           NumFilters;         // number of filter records
    short           NumCoeffs;          // number of coefficient records
    int16_t         Stmp;               // short temporary variable
    filter          *Filter;            // Filter[F] -- array of filters
    NewSensorResRec SensorRecord;       // sensor information
    SensorCoefResRec CoeffRecord;       // ugly balance coefficient info records

    // split path into dspath and setname, save in header
    if(ParsePath(DsName, DsPath) == -1)
        cleanup("'ParsePath()' failed");
    Header->DsPath = (char *)malloc(strlen(DsPath) + 1);
    strcpy(Header->DsPath, DsPath);
    Header->SetName = (char *)malloc(strlen(DsName) + 1);
    strcpy(Header->SetName, DsName);

    // check for access!
    sprintf(PathName, "%s/%s.ds", Header->DsPath, Header->SetName);
    if(access(PathName, F_OK) != 0)
        cleanup("dataset %s not found!", PathName);

    // 1st, read resource file
    sprintf(PathName, "%s/%s.ds/%s.res4", Header->DsPath, Header->SetName, Header->SetName);
    if((fin = fopen(PathName, "r")) == NULL)
        cleanup("'.res4' file open failed");

    // read & verify header label
    if(fread((void *)&HeaderLabel, 8, 1, fin) != 1)
        cleanup("'HeaderLabel' read failed");
#if 0
    if(!strncmp(HeaderLabel, "MEG42", 5)) {
        MEGVersion = 42;
    } else if(!strncmp(HeaderLabel, "MEG41", 5)) {
        MEGVersion = 41;
    } else if(!strncmp(HeaderLabel, "MEG40", 5)) {
        MEGVersion = 40;
    } else if(!strncmp(HeaderLabel, "MEG4RES", 7)) {
        MEGVersion = 40;
    } else {
        cleanup("unknown file ID");q
    }
#endif

    // read general resource & fill 'HeaderInfo' structure
    if(fread((void *)&GeneralResource, 1832, 1, fin) != 1)
        cleanup("'GeneralResource' read failed");
    if(fread((void *)&Padding, 4, 1, fin) != 1)
        cleanup("'Padding read failed");
    NumBytes = bswap_32(GeneralResource.nfSetUp.size);
    if((run_descriptor = (char *)malloc(NumBytes)) == NULL)
        allocfailed("run_descriptor");
    if(fread((void *)run_descriptor, NumBytes, 1, fin) != 1)
        cleanup("'run_descriptor' read failed");
    if(fread((void *)&Stmp, sizeof(short int), 1, fin) != 1)
        cleanup("'NumFilters' read failed");
    NumFilters = bswap_16(Stmp);
    if((Filter = (filter *)malloc((24 * NumFilters))) == NULL)
        allocfailed("Filter");
    for(i=0; i<NumFilters; i++) {
        if(fread((void *)&Dtmp, sizeof(double), 1, fin) != 1)
            cleanup("'freq' read failed");
        Filter[i].freq = Dflip(&Dtmp);
        if(fread((void *)&Itmp, sizeof(int), 1, fin) != 1)
            cleanup("'fClass' read failed");
        Filter[i].fClass = (classType)bswap_32(Itmp);
        if(fread((void *)&Itmp, sizeof(int), 1, fin) != 1)
            cleanup("'fType' read failed");
        Filter[i].fType = (filtType)bswap_32(Itmp);
        if(fread((void *)&Stmp, sizeof(short int), 1, fin) != 1)
            cleanup("'numParam' read failed");
        Filter[i].numParam = bswap_16(Stmp);
        if((Filter[i].params = (double *)malloc( (8 * Filter[i].numParam))) == NULL)
            allocfailed("Filter[i].params");
        for(j=0; j<Filter[i].numParam; j++) {
            if(fread((void *)&Dtmp, sizeof(double), 1, fin) != 1)
                cleanup("'params' read failed");
            Filter[i].params[j] = Dflip(&Dtmp);
        }
    }

    // create 'Header' structure
    Header->Version = 1;
    M = (int)bswap_16(GeneralResource.gSetUp.no_channels);
    Header->NumChannels = Header->NumAcquired = M;
    E = (int)bswap_16(GeneralResource.gSetUp.no_trials);
    Header->NumEpochs = E;
    T = bswap_32(GeneralResource.gSetUp.no_samples);
    Header->MaxSamples = T;
    Header->SampleRate = Dflip(&GeneralResource.gSetUp.sample_rate);
    Header->PowerFreq = 60.;    // eh?
    Header->PreTrigSamples = (int)bswap_32(GeneralResource.gSetUp.preTrigPts);
    Header->Avg = bswap_16(GeneralResource.no_trials_avgd);
    Header->DataFormat = 0;
    Header->spillt = (MAX_MEG4 >> 2) / (M * T);     // spillt has the theoretical maximum epoch count in a .meg4 file

    // allocate & read channel information
    if((Channel = (ChannelInfo *)malloc(M * sizeof(ChannelInfo))) == NULL)
        allocfailed("Channel");
    if((Gain = (double *)malloc(M * sizeof(double))) == NULL)
        allocfailed("Gain");
    *ChannelOut = Channel;      // return pointer to channel info array

    // read channel names
    for(m=0; m<M; m++) {
        if(fread((void *)ChannelLabel, 32, 1, fin) != 1)
            cleanup("'ChannelLabel' read failed");
        l = strlen(ChannelLabel);
        memcpy(Channel[m].ChannelName, ChannelLabel, l+1);
    }

    // read channel resources
    for(m=NumSens=NumRefs=NumEEGs=0; m<M; m++) {
        if(fread((void *)&SensorRecord, sizeof(NewSensorResRec), 1, fin) != 1)
            cleanup("'SensorRecord' read failed");
        switch(bswap_16(SensorRecord.sensorTypeIndex)) {
            case eMEGReference:     // MEG magnetometer references
                Channel[m].ChannelType = TYPE_REF_MAG;
                Channel[m].Units = UNIT_TESLA;
                NumRefs++;          // count reference channels
                break;
            case eMEGReference1:    // MEG gradiometer references
            case eMEGReference2:
            case eMEGReference3:
                Channel[m].ChannelType = TYPE_REF_GRAD;
                Channel[m].Units = UNIT_TESLA;
                NumRefs++;          // count reference channels
                break;
            case eMEGSensor:        // MEG gradiometer or magnetometer sensors
            case eMEGSensor1:
            case eMEGSensor2:
            case eMEGSensor3:
                Channel[m].ChannelType = TYPE_MEG;
                Channel[m].Units = UNIT_TESLA;
                NumSens++;          // count sensor channels
                break;
            case eEEGRef:           // EEG reference channels
                Channel[m].ChannelType = TYPE_OTHER;
                Channel[m].Units = UNIT_VOLTS;
                break;
            case eEEGSensor:        // EEG sensor channels
                Channel[m].ChannelType = TYPE_EEG;
                Channel[m].Units = UNIT_VOLTS;
                NumEEGs++;
                break;
            case eADCRef:           // ADC reference channels
                Channel[m].ChannelType = TYPE_EXTERNAL;
                Channel[m].Units = UNIT_BITS;
                break;
            case eStimRef:          // Stimulus/trigger channel
                Channel[m].ChannelType = TYPE_TRIGGER;
                Channel[m].Units = UNIT_BITS;
                //TriggerChannel = m;   // save trigger channel index
                break;
            default:
                Channel[m].ChannelType = TYPE_OTHER;
                Channel[m].Units = UNIT_OTHER;
                break;
        }
        Channel[m].Flag = TRUE;
        Channel[m].ChannelNumber = m;
        qGain = Dflip(&SensorRecord.qGain);
        ioGain = Dflip(&SensorRecord.ioGain);
        properGain = Dflip(&SensorRecord.properGain);
        Channel[m].Scale = 1. / (qGain * ioGain * properGain);
        Channel[m].AnalogHpFreq = 0;
        Channel[m].AnalogLpFreq = Header->SampleRate / 4.;
        for(i=0; i<NumFilters; i++)
            if(Filter[i].fClass == BUTTERWORTH && Filter[i].fType == LOWPASS) {
                Channel[m].AnalogLpFreq = Filter[i].freq;
            break;
        }
        for(i=0; i<NumFilters; i++)
            if(Filter[i].fClass == BUTTERWORTH && Filter[i].fType == HIGHPASS) {
                Channel[m].AnalogHpFreq = Filter[i].freq;
            break;
        }
        Gain[m] = properGain;       // save 'properGain' for correcting balance
        Polarity = Gain[m] / fabs(Gain[m]);
        if(Channel[m].ChannelType == TYPE_MEG || Channel[m].ChannelType == TYPE_REF_MAG || Channel[m].ChannelType == TYPE_REF_GRAD) {

            // for magnetic sensors, find out how many coils per sensor & allocate coil structures
            Channel[m].Geom.MEGSensor.NumCoils = bswap_16(SensorRecord.numCoils);
            if((Channel[m].Geom.MEGSensor.Coil = (COIL *) malloc( (Channel[m].Geom.MEGSensor.NumCoils * sizeof(COIL)))) == NULL)
                allocfailed("Coil");

            // first, fill in COIL parameters from CTF format
            for(c=0; c<Channel[m].Geom.MEGSensor.NumCoils; c++) {
            if(Xformed == TRUE) {
                for(v=X_; v<=Z_; v++) {
                    Channel[m].Geom.MEGSensor.Coil[c].origin.p[v] = 0.01 * Dflip(&SensorRecord.HdcoilTbl[c].position.point[v]);
                    Channel[m].Geom.MEGSensor.Coil[c].origin.v[v] = Polarity * Dflip(&SensorRecord.HdcoilTbl[c].orient.point[v]);
                }
            } else {
                for(v=X_; v<=Z_; v++) {
                    Channel[m].Geom.MEGSensor.Coil[c].origin.p[v] = 0.01 * Dflip(&SensorRecord.coilTbl[c].position.point[v]);
                    Channel[m].Geom.MEGSensor.Coil[c].origin.v[v] = Polarity * Dflip(&SensorRecord.coilTbl[c].orient.point[v]);
                }
            }
            Channel[m].Geom.MEGSensor.Coil[c].radius = 0.01 * sqrt(Dflip(&SensorRecord.coilTbl[c].area) / M_PI);
            Channel[m].Geom.MEGSensor.Coil[c].SenseTurns = bswap_16(SensorRecord.coilTbl[c].numturns);
            for(i=0; i<INTPTS; i++)
                for(v=X_; v<=Z_; v++) {
                    Channel[m].Geom.MEGSensor.Coil[c].B[i].p[v] = 0.;
                    Channel[m].Geom.MEGSensor.Coil[c].B[i].v[v] = 0.;
                }
                Channel[m].Geom.MEGSensor.Coil[c].Evaluated = FALSE;
            }

            // next, fill in SENSOR structure
            for(v=X_; v<=Z_; v++) {
                Channel[m].Geom.MEGSensor.origin.p[v] = Channel[m].Geom.MEGSensor.Coil[0].origin.p[v];
                Channel[m].Geom.MEGSensor.origin.v[v] = Channel[m].Geom.MEGSensor.Coil[0].origin.v[v];
                Channel[m].Geom.MEGSensor.LocalSphere[v] = 0.;
            }

            // set up balance coefficients on primary sensors (as long as refs come before primaries!)
            if (Channel[m].ChannelType == TYPE_MEG) {
                Channel[m].Balance = new_arrayE(BalInfo, NumRefs, "Balance[]");
                for (n=0; n<NumRefs; n++) {
                    Channel[m].Balance[n].Weight = 0.;
                    Channel[m].Balance[n].RefChan = -1;
                }
                Channel[m].NumBal = 0;
                Channel[m].BalOrder = bswap_16(SensorRecord.grad_order_no);
            } else {
                Channel[m].Balance = NULL;      // reference channels have no balance info
                Channel[m].NumBal = 0;
                Channel[m].BalOrder = -1;
            }
        } else {
            Channel[m].Geom.NoGeom = TRUE;
            Channel[m].Balance = NULL;
            Channel[m].NumBal = 0;
            Channel[m].BalOrder = -1;
        }
    }

    // attempt to read optional .eeg file
    sprintf(PathName, "%s/%s.ds/%s.eeg", Header->DsPath, Header->SetName, Header->SetName);
    if((feeg = fopen(PathName, "r")) != NULL) {

        // allocate space for PeIndex
        Header->NumEEG = NumEEGs;
        if((Header->PeIndex = (int *)malloc(Header->NumEEG * sizeof(int))) == NULL)
            allocfailed("PeIndex[]");

        // read all lines of .eeg file
        for(i=0; i<Header->NumEEG; i++) {
            if(fscanf(feeg, "%d%s%lf%lf%lf", &n, ChannelLabel, &Vc[X_], &Vc[Y_], &Vc[Z_]) != 5)
                cleanup("failed to read line from .eeg file");
            m = n - 1;
            if(m >= Header->NumAcquired)
                cleanup("eeg channel index exceeds total channels in 'res4' file");
            strncpy(Channel[m].ChannelName, ChannelLabel, sizeof(Channel[m].ChannelName));
            for(v=X_; v<=Z_; v++)
                Channel[m].Geom.EEGSensor.position[v] = 0.01 * Vc[v];
            Channel[m].Geom.EEGSensor.gain = 1.;
            Channel[m].Geom.EEGSensor.RefNumber = m;
            Channel[m].Balance = NULL;
            Channel[m].NumBal = 0;
            Channel[m].BalOrder = -1;
            Header->PeIndex[i] = n;
        }
        fclose(feeg);
    }

    // read balance coefficients
    if (fread((void *)&Stmp, sizeof(int16_t), 1, fin) != 1)
        cleanup("'NumCoeffs' read failed");
    NumCoeffs = bswap_16(Stmp);
    for (i=0; i<NumCoeffs; i++) {       // there should be at least one of these per channel/order

        // read one sensor record
        if (fread((void *)&CoeffRecord, 1992, 1, fin) != 1)
            cleanup("'CoeffRecord' read failed");

        // find SQUID sensor channel index m
        l = strlen(CoeffRecord.sensorName);
        for (m=0, flag=FALSE; m<M; m++) // convert sensor label to sensor channel index
            if (!strncmp(Channel[m].ChannelName, CoeffRecord.sensorName, l)) {
                flag = TRUE;
                break;
            }
        if (!flag)                      // make certain that we found the channel index
            cleanup("channel index not found");

        // determine the balance order to be used in forward solutions
        order = -1;
        switch (Channel[m].BalOrder) {
            case -1:
                order = -1;
                break;
            case 0:
                order = NOGRAD;
                break;
            case 1:
                order = G1BR;
                break;
            case 2:
                order = G2BR;
                break;
            case 3:
                order = G3BR;
                break;
            case 10:
                order = G0AR;
                break;
            case 11:
                order = G1AR;
                break;
            case 12:
                order = G2AR;
                break;
            case 13:
                order = G3AR;
                break;
            default:
                cleanup("unknown balance order");
        }

        // now, determine if the balance order for this record matches the stored order
        if (order != -1 && bswap_32(CoeffRecord.coefType) == order) {

            // save number of coefficients, & for each coefficient j of sensor m...
            Channel[m].NumBal = bswap_16(CoeffRecord.coefRec.num_of_coefs);
            for (j=0; j<Channel[m].NumBal; j++) {

                // find SQUID reference channel n
                l = strlen(CoeffRecord.coefRec.sensor_list[j]);
                for (n=0, flag=FALSE; n < M; n++)   // convert reference label to channel index
                    if (!strncmp(Channel[n].ChannelName, CoeffRecord.coefRec.sensor_list[j], l)) {
                        flag = TRUE;
                        break;
                    }
                if (flag) {                         // make certain that we found the channel index
                    // make certain that this is a reference
                    if (Channel[n].ChannelType != TYPE_REF_MAG && Channel[n].ChannelType != TYPE_REF_GRAD)
                        cleanup("indexed channel was not a reference");
                    Channel[m].Balance[j].Weight = Dflip(&CoeffRecord.coefRec.coefs_list[j]) * Gain[n] / Gain[m];
                    Channel[m].Balance[j].RefChan = n;
                } else {
                    break;
                }
            }
        }
    }

    // allocate & read epoch information
    if((Epoch = (EpochInfo *)malloc(E * sizeof(EpochInfo))) == NULL)
        allocfailed("Epoch");
    *EpochOut = Epoch;          // return pointer to epoch info array
    for(e=0; e<E; e++) {
        Epoch[e].NumSamples = T;
        Epoch[e].Flag = TRUE;   // for now...
//      Itmp = bswap_32(GeneralResource.gSetUp.preTrigPts);
//      Epoch[e].StimulusTime = (float)Itmp / Header->SampleRate;   // 'getSamplePreTrigger()' returns sample number!
    }
    fclose(fin);

    // open 'ClassFile.cls', & check if there are 'BAD' epochs
    sprintf(PathName, "%s/%s.ds/ClassFile.cls", Header->DsPath, Header->SetName);
    if((fp = fopen(PathName, "r")) != NULL) {

        do {                    // synch with & read number of classes
            if(fgets(StringField, 128, fp) == NULL)
            goto nobad;
        } while(strncmp(StringField, "NUMBER OF CLASSES:", 17));
        if(fscanf(fp, "%d", &NumClasses) == EOF)
            goto nobad;

        // scan each class to find 'BAD' class name
        for(i=0; i<NumClasses; i++) {
            do {
                if(fgets(StringField, 128, fp) == NULL)
                    goto nobad;
            } while(strncmp(StringField, "NAME:", 5));
            if(fgets(StringField, 128, fp) == NULL)
                goto nobad;
            if(!strncmp(StringField, "BAD", 3)) {

                do {            // search for number of trials in 'BAD' class
                    if(fgets(StringField, 64, fp) == NULL)
                        goto nobad;
                } while(strncmp(StringField, "NUMBER OF TRIALS:", 17));
                if(fscanf(fp, "%d", &NumBad) == EOF)
                    goto nobad;

                do {            // synch up to line before marker list
                    if(fgets(StringField, 64, fp) == NULL)
                    goto nobad;
                } while(strncmp(StringField, "TRIAL NUMBER", 12));

                // read list of epochs marked 'BAD'
                for(i=0; i<NumBad; i++) {
                    if(fscanf(fp, "%d", &j) == EOF)
                        goto nobad;
                    Epoch[j].Flag = FALSE;
                }
            }
        }
        nobad:
        fclose(fp);
    }

    // open 'BadChannels' & read names of bad channels
    sprintf(PathName, "%s/%s.ds/BadChannels", Header->DsPath, Header->SetName);
    if((fp = fopen(PathName, "r")) != NULL) {
        while(fgets(StringField, 128, fp) != NULL) {
            NumBytes = strlen(StringField) - 1;
            for(m=0; m<M; m++)
                if(!strncmp(StringField, Channel[m].ChannelName, NumBytes)) {
                    Channel[m].Flag = FALSE;
                    break;
                }
        }
        fclose(fp);
    }

    // save number of primary & reference sensors with indices
    Header->NumPri = NumSens;
    Header->NumRef = NumRefs;
    Header->NumSquid = NumRefs + NumSens;
    if((Header->PsIndex = (int *)malloc(Header->NumPri * sizeof(int))) == NULL)
        allocfailed("PsIndex");
    if((Header->RsIndex = (int *)malloc(Header->NumRef * sizeof(int))) == NULL)
        allocfailed("RsIndex");
    if((Header->SqIndex = (int *)malloc(Header->NumSquid * sizeof(int))) == NULL)
        allocfailed("SqIndex");

    // get primary SQUID indices
    for(m=i=0; m<M; m++)
        if(Channel[m].ChannelType == TYPE_MEG) {
            if(Channel[m].Flag == TRUE) {
                Header->PsIndex[i] = m;
                i++;
            } else {
                Header->NumPri--;
            }
        }

    // get reference SQUID indices
    for(m=i=0; m<M; m++)
        if(Channel[m].ChannelType == TYPE_REF_MAG || Channel[m].ChannelType == TYPE_REF_GRAD) {
            if(Channel[m].Flag == TRUE) {
                Header->RsIndex[i] = m;
                i++;
            } else {
                Header->NumRef--;
            }
        }

    // get all SQUID indices
    for(m=i=0; m<M; m++)
        if(Channel[m].ChannelType == TYPE_MEG || Channel[m].ChannelType == TYPE_REF_MAG || Channel[m].ChannelType == TYPE_REF_GRAD) {
            if(Channel[m].Flag == TRUE) {
                Header->SqIndex[i] = m;
                i++;
            } else {
                Header->NumSquid--;
            }
        }

    // if req'd, read 'bad.samples' & create array of flags
    if(BadOut != NULL) {

        // allocate 'Bad' samples
        if((Bad = (unsigned char **)malloc(E * sizeof(unsigned char *))) == NULL)
            allocfailed("Bad[]");
        for(e=0; e<E; e++)
            if((Bad[e] = (unsigned char *)malloc(T * sizeof(unsigned char))) == NULL)
                allocfailed("Bad[][]");

        // initialize bad samples to agree with BAD epochs
        for(e=0; e<E; e++) {
            flag = (Epoch[e].Flag == TRUE) ? GOOD : BAD;
            for(t=0; t<T; t++)
                Bad[e][t] = (unsigned char)flag;
        }

        // if present, open 'bad.segments file, set bad samples, & return pointer
        sprintf(PathName, "%s/%s.ds/bad.segments", Header->DsPath, Header->SetName);
        if((fp = fopen(PathName, "r")) != NULL) {       // bad.segments exists!
            while(fscanf(fp, "%d%lf%lf", &BadEpoch, &BadStart, &BadEnd) == 3) {
                e = BadEpoch - 1;                       // 'cause the epoch numbers are base 0
                if(e < 0 || e > E-1)                    // check for error in file
                    cleanup("epoch number out of range in 'bad.segments'");
                bs = (int)rint(BadStart * Header->SampleRate);
                if(bs < 0)
                    bs = 0;
                be = (int)rint(BadEnd * Header->SampleRate);
                if(be >= T)
                    be = T - 1;
                for(t=bs; t<= be; t++)
                    Bad[e][t] = BAD;
            }
            fclose(fp);
        }
        *BadOut = Bad;          // remember to return pointer to 'BadSample' -- that's all!
    }
}
