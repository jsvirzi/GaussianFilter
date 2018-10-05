#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

#include <TApplication.h>
#include <TGraph.h>
#include <TCanvas.h>

int main(int argc, char **argv) {
    size_t n = 1024;
    double dt = 0.001;
    unsigned int xFilterLength = 61, yFilterLength = 2;
    double T = xFilterLength * dt;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-n") == 0) {
            n = atoi(argv[++i]);
        }
    }

    TApplication the_app("theApp", &argc, argv);

    double *xcra = new double [xFilterLength];
    double *xcrb = new double [xFilterLength];
    double *xci = new double [xFilterLength];
    double *ycr = new double [yFilterLength];
    double *yci = new double [yFilterLength];

    double omega0 = 2.0 * M_PI * 100.0;
    double a = 1.0 / omega0; /* for proper normalization */
    double sigma = 2.0 * M_PI * 5.0; /* 5 Hz half-width */
    unsigned int i;
    for (i = 0; i < xFilterLength / 2; ++i) {
        double t = i * dt - 0.5 * T;
        double b = sigma * sigma * t * t;
        xcra[i] = a * sin(omega0 * t) / t;
        xcrb[i] = xcra[i] * exp(-1.0 * b);
        xcra[xFilterLength - 1 - i] = xcra[i];
        xcrb[xFilterLength - 1 - i] = xcrb[i];
    }
    xcra[i] = 1.0;
    xcrb[i] = 1.0;

    unsigned int nSeconds = 5;
    unsigned int samplingFrequency = 1000;
    double fs = samplingFrequency;
    double tStart = 0.0, tFinal = nSeconds;
    unsigned int arraySize = nSeconds * samplingFrequency;
    double *si = new double [arraySize];
    double *sa = new double [arraySize];
    double *sb = new double [arraySize];
    double *tt = new double [arraySize];
    for (unsigned int i = 0; i < arraySize; ++i) {
        double t = i * dt;
        tt[i] = t;
    }

    unsigned int nFrequencies = 200;
    double *fa = new double [nFrequencies];
    double *fb = new double [nFrequencies];
    double *fx = new double [nFrequencies];
    for (int j = 0; j < nFrequencies; ++j) {
        fx[j] = j;
    }

    for (int j = 0; j < nFrequencies; ++j) {
        double omega = 2.0 * M_PI * fx[j];

        for (unsigned int i = 0; i < arraySize; ++i) {
            double t = tt[i];
            si[i] = sin(omega * t);
            sa[i] = 0.0;
            sb[i] = 0.0;
        }

        for (unsigned int i = xFilterLength / 2; i < arraySize; ++i) {
            double acca = 0.0, accb = 0.0;
            for (unsigned int j = 0; j < xFilterLength; ++j) {
                acca = acca + xcra[j] * si[i + xFilterLength / 2 - j];
                accb = accb + xcrb[j] * si[i + xFilterLength / 2 - j];
            }
            sa[i] = acca;
            sb[i] = accb;
        }

        double ea = 0.0, eb = 0.0;
        for (unsigned int i = xFilterLength / 2; i < arraySize; ++i) {
            ea = ea + sa[i] * sa[i];
            eb = eb + sb[i] * sb[i];
        }
        fa[j] = ea;
        fb[j] = eb;
    }

    TGraph *g, *ga, *gb;

    int nDisplay = nSeconds * 25;
    double *tDisplay = new double [nDisplay];
    double *yDisplay = new double [nDisplay];

    tDisplay[0] = 0.0;
    yDisplay[0] = 1.0;
    for (int i = 1; i < nDisplay; ++i) {
        int index = 10 * i;
        tDisplay[i] = tt[index];
        yDisplay[i] = sa[index];
    }

    TGraph *graphXCRa = new TGraph(xFilterLength, tt, xcra);
    TGraph *graphXCRb = new TGraph(xFilterLength, tt, xcrb);
    TGraph *graphYa = new TGraph(arraySize, tt, sa);
    TGraph *graphYb = new TGraph(arraySize, tt, sb);
    TGraph *graphEa = new TGraph(nFrequencies, fx, fa);
    TGraph *graphEb = new TGraph(nFrequencies, fx, fb);
//    TGraph *graphXCR = new TGraph(nDisplay, tDisplay, xcra);

    ga = graphEa;
    gb = graphEb;

    ga->SetLineColor(kBlack);
    ga->SetLineWidth(3.0);
    gb->SetLineColor(kRed);
    gb->SetLineWidth(3.0);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 3000, 600);
    ga->Draw("AL");
    gb->Draw("L");
    canvas->Update();
    canvas->Draw();
    canvas->WaitPrimitive();

    return 0;
}
