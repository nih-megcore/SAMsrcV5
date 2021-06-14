#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lfulib.h>
#include <samlib.h>

int CoeffFlag;          // @@@ delete me
#define NUM_INT     3           // number of integration points per coil (<= INTPTS)

// This is a (not thread-safe!) convenience function for programs that only use one sensor frame.
// SAMcoreg needs to call ECDIntPnt/ECDIntPntFree every time the sensor moves.

void ECDLeadField(
    double      Org[3],         // ECD coordinate
    gsl_matrix  *L,             // L - Mx3 output ECD primary sensor lead field matrix
    HeaderInfo  *Header,        // MEG dataset header information
    ChannelInfo *Channel,       // Channel[M] -- MEG channel information
    HULL        *Hull           // cortical hull model
) {
    static COEFFS *c = NULL;    // structure to store detection points and coefficients (per sensor frame)

    if (c == NULL) {
        c = new(COEFFS);
        ECDIntPnt(Header, Channel, Hull, c);
    }
    ECDLeadField0(Org, L, Header, Channel, Hull, c);
}

// Compute detection points and coefficents, when the sensor frame is first established, or if it changes.
// Remember to call ECDIntPntFree() before calling this again!

void ECDIntPnt(
    HeaderInfo  *Header,        // MEG dataset header information
    ChannelInfo *Channel,       // Channel[M] -- MEG channel information
    HULL        *Hull,          // cortical hull model
    COEFFS      *coeff          // space to return detection point and coefficient info
) {
    int         i, m, n, c, v;
    int         S, C;           // number of integration points and coils
    int         BASIS;          // number of basis functions
    double      x[3];           // x' vector
    double      y[3];           // y' vector
    double      Radius;         // integration point radius
    double      angle;          // integration point in-plane angle
    double      XDisp;          // x-displacement
    double      YDisp;          // y-displacement
    FIELD       *detector;      // detector[S] -- integration points
    double      **coeffs;       // coeffs[BASIS][S] -- correction coefficients

    // count integration points & coils

    for (m=S=C=0; m<Header->NumChannels; m++)
        if ((Channel[m].ChannelType == TYPE_REF_MAG ||
             Channel[m].ChannelType == TYPE_REF_GRAD ||
             Channel[m].ChannelType == TYPE_MEG) && Channel[m].Flag)
            for (c=0; c<Channel[m].Geom.MEGSensor.NumCoils; c++) {
                C += 1;
                S += NUM_INT;
            }

    BASIS = ((ORDER+1) * (ORDER+1) - 1);
    coeff->C = C;
    coeff->S = S;
    coeff->detector = detector = new_arrayE(FIELD, S, "detector[]");
    coeff->coeffs = coeffs = new_matrixE(BASIS, S, "coeffs[][]");

    for (m=0; m<Header->NumChannels; m++) {
        if ((Channel[m].ChannelType == TYPE_REF_MAG ||
             Channel[m].ChannelType == TYPE_REF_GRAD ||
             Channel[m].ChannelType == TYPE_MEG) && Channel[m].Flag) {

            // evaluate all coil integration points
            for (c=0; c<Channel[m].Geom.MEGSensor.NumCoils; c++) {

                // set up three integration points in plane of coil
                OrthoPlane(Channel[m].Geom.MEGSensor.Coil[c].origin.v, x, y);   // solve for x' & y' defining plane of coil
                Channel[m].Geom.MEGSensor.Coil[c].NumInt = NUM_INT;
                Radius = M_SQRT1_2 * Channel[m].Geom.MEGSensor.Coil[c].radius;
                for (i=0; i<NUM_INT; i++) {
                    angle = (double)i * 2. * M_PI / (double)NUM_INT;
                    XDisp = Radius * cos(angle);
                    YDisp = Radius * sin(angle);
                    for (v=X_; v<=Z_; v++)
                        Channel[m].Geom.MEGSensor.Coil[c].B[i].p[v] =
                            Channel[m].Geom.MEGSensor.Coil[c].origin.p[v] + x[v] * XDisp + y[v] * YDisp;
                }
                Channel[m].Geom.MEGSensor.Coil[c].Evaluated = TRUE;
            }
        }
    }

    // copy detection points to array
    for (m=n=0; m<Header->NumChannels; m++) {
        if ((Channel[m].ChannelType == TYPE_REF_MAG ||
             Channel[m].ChannelType == TYPE_REF_GRAD ||
             Channel[m].ChannelType == TYPE_MEG) && Channel[m].Flag) {
            for (c=0; c<Channel[m].Geom.MEGSensor.NumCoils; c++) {
                for (i=0; i<NUM_INT; i++) {
                    for (v=X_; v<=Z_; v++) {
                        detector[n].p[v] = Channel[m].Geom.MEGSensor.Coil[c].B[i].p[v];
                        detector[n].v[v] = Channel[m].Geom.MEGSensor.origin.v[v];
                    }
                    n++;
                }
            }
        }
    }

    // compute coefficients
    GetCoeffs(detector, S, Hull, coeffs);
}

void ECDIntPntFree(
    COEFFS      *c          // detection point info
) {
    free(c->detector);
    free_matrix(c->coeffs);
}

void ECDLeadField0(
    double      Org[3],                 // ECD coordinate
    gsl_matrix  *L,                     // L - Mx3 output ECD primary sensor lead field matrix
    HeaderInfo  *Header,                // MEG dataset header information
    ChannelInfo *Channel,               // Channel[M] -- MEG channel information
    HULL        *Hull,                  // cortical hull model
    COEFFS      *coeff                  // detection point info
) {
    FIELD       *detector;              // detector[S] -- integration points
    double      **coeffs;               // coeffs[BASIS][S] -- correction coefficients
    double      **gradbasis;            // gradbasis[BASIS+1][3] -- gradient of basis functions
    double      **Bi;                   // Bi[S][3] -- corrected forward solution for detection points
    double      Bd[3][3];
    double      Vd[3];                  // dipole location
    double      Vs[3];                  // detector location
    double      **fieldori;             // fieldori[S][3]
    double      **correction;           // correction[S][3]
    double      **coil;                 // coil[C][3]
    double      **Bp;                   // Bp[M][3] -- primary sensor lead field
    double      **Br;                   // Br[R][3] -- reference sensor lead field
    int         c;                      // coil index
    int         cc;                     // alternate coil index
    int         g;
    int         i, j;
    int         m;                      // sensor index
    int         found;                  // found reference channel
    int         v, w;                   // vector indices
    int         BASIS;                  // number of basis functions
    int         C;                      // number of coils
    int         S;                      // number of integration points

    // extract passed coefficient data

    C = coeff->C;
    S = coeff->S;
    detector = coeff->detector;
    coeffs = coeff->coeffs;

    // allocate memory

    BASIS = ((ORDER+1) * (ORDER+1) - 1);
    gradbasis = new_matrixE(BASIS + 1, 3, "gradbasis[][]");
    fieldori = new_matrixE(S, 3, "fieldori[][]");
    correction = new_matrixE(S, 3, "correction[][]");
    Bi = new_matrixE(S, 3, "Bi[][]");
    coil = new_matrixE(C, 3, "coil[][]");
    Bp = new_matrixE(Header->NumPri, 3, "Bp[][]");
    Br = NULL;
    if (Header->NumRef > 0) {
        Br = new_matrixE(Header->NumRef, 3, "Br[][]");
    }

    // translate dipole source location to volume origin & compute correction for this dipole position
    for (v=X_; v<=Z_; v++)
        Vd[v] = Org[v] - Hull->Vo[v];
    GetBasis(Vd, Hull->scale, gradbasis);

    // for all detection points
    for (i=0; i<S; i++) {

        // translate detector location to conducting volume origin
        for (v=X_; v<=Z_; v++)
            Vs[v] = detector[i].p[v] - Hull->Vo[v];

        // solve field for spherical volume conductor (origin of volume conductor)
        xbd(Vd, Vs, Bd);

        // for each cardinal moment (v)
        for (v=X_; v<=Z_; v++) {

            // field is dot product of field Bd with detector normal vector (w)
            for (w=X_, fieldori[i][v]=0.; w<=Z_; w++)
                fieldori[i][v] += Bd[v][w] * detector[i].v[w];

            // final result
            for (j=0, correction[i][v]=0.; j<BASIS; j++)
                correction[i][v] -= gradbasis[j+1][v] * coeffs[j][i];       // note: gradbasis 1...BASIS
            Bi[i][v] = fieldori[i][v] + correction[i][v];
        }
    }

    // integrate field within each coil
    for (i=0, c=0; i<S; i+=NUM_INT, c++)
        for (v=X_; v<=Z_; v++)
            for (j=0, coil[c][v]=0.; j<NUM_INT; j++)
                coil[c][v] += Bi[i+j][v] / (double)NUM_INT;

    // compute field for each sensor
    for (m=c=i=j=0; m<Header->NumChannels; m++)
        if ((Channel[m].ChannelType == TYPE_REF_MAG ||
             Channel[m].ChannelType == TYPE_REF_GRAD) && Channel[m].Flag) {
            for (v=X_; v<=Z_; v++) {
                for (cc=0, Br[i][v]=0.; cc<Channel[m].Geom.MEGSensor.NumCoils; cc++)
                    Br[i][v] += coil[c][v] * (double)Channel[m].Geom.MEGSensor.Coil[cc].SenseTurns;
                Br[i][v] /= (double)abs(Channel[m].Geom.MEGSensor.Coil[0].SenseTurns);
            }
            c += Channel[m].Geom.MEGSensor.NumCoils;
            i++;
        } else if (Channel[m].ChannelType == TYPE_MEG && Channel[m].Flag) {
            for (v=X_; v<=Z_; v++) {
                for (cc=0, Bp[j][v]=0.; cc<Channel[m].Geom.MEGSensor.NumCoils; cc++)
                    Bp[j][v] += coil[c][v] * (double)Channel[m].Geom.MEGSensor.Coil[cc].SenseTurns;
                Bp[j][v] /= (double)abs(Channel[m].Geom.MEGSensor.Coil[0].SenseTurns);
            }
            c += Channel[m].Geom.MEGSensor.NumCoils;
            j++;
        }

    // if req'd, apply balance coefficients & output primary sensor lead field
    if (Header->NumRef > 0) {   // only if there are reference sensors
        for (m=0; m<Header->NumPri; m++) {
            j = Header->PsIndex[m];
            for (g=0; g<Channel[j].NumBal; g++) {
                for (i=0, found=FALSE; i<Header->NumRef; i++)
                    if (Channel[j].Balance[g].RefChan == Header->RsIndex[i]) {
                        found = TRUE;
                        break;
                    }
                if (!found)
                    cleanup("couldn't find reference channel index");
                for (v=X_; v<=Z_; v++)
                    Bp[m][v] -= Channel[j].Balance[g].Weight * Br[i][v];
            }
        }
    }
    for (m=0; m<Header->NumPri; m++)
        for (v=X_; v<=Z_; v++)
            gsl_matrix_set(L, m, v, Bp[m][v]);

    free_matrix(gradbasis);
    free_matrix(fieldori);
    free_matrix(correction);
    free_matrix(Bi);
    free_matrix(coil);
    free_matrix(Bp);
    if (Br) {
        free_matrix(Br);
    }
}
