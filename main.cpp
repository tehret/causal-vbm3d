#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>
#include <algorithm>

#include "vbm3d.h"
#include "vpp/vpp.h"
#include "Utilities/Utilities.h"
#include "Utilities/cmd_option.h"
#include "lib_transforms.h"

#define YUV       0
#define YCBCR     1
#define OPP       2
#define RGB       3
#define DCT       4
#define BIOR      5
#define HADAMARD  6
#define HAAR      7
#define NONE      8

using namespace std;

void initializeParameters_1(
	Parameters& prms
,	const int k
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float tau
,	const float lambda3D
,	const unsigned T_2D
,	const unsigned T_3D
,	const float sigma
){
	if(k < 0)
		prms.k = 8;
	else
		prms.k = k;

	if(Nf < 0)
		prms.Nf = 4;
	else
		prms.Nf = Nf;
	if(Ns < 0)
		prms.Ns = 7;
	else
		prms.Ns = Ns;
	if(Npr < 0)
		prms.Npr = 5;
	else
		prms.Npr = Npr;
	if(Nb < 0)
		prms.Nb = 2;
	else
		prms.Nb = Nb;

	if(d < 0)
		prms.d = (7*7)/(255.*255.);
	else
		prms.d = (d*d)/(255.*255.);

	if(p < 0)
		prms.p = 6;
	else
		prms.p = p;

	if(N < 0)
		prms.N = 8;
	else
		prms.N = N;

	if(tau < 0)
		prms.tau = (sigma > 30) ? 4500/(255.*255.) : 3000/(255.*255.);
	else
		prms.tau = tau;

	if(lambda3D < 0)
		prms.lambda3D = 2.7f;
	else
		prms.lambda3D = lambda3D;

	if(T_2D == NONE)
		prms.T_2D = BIOR;
	else
		prms.T_2D = T_2D;

	if(T_3D == NONE)
		prms.T_3D = HAAR;
	else
		prms.T_3D = T_3D;
}

void initializeParameters_2(
	Parameters& prms
,	const int k
,	const int Nf
,	const int Ns
,	const int Npr
,	const int Nb
,	const int p
,	const int N
,	const int d
,	const float tau
,	const unsigned T_2D
,	const unsigned T_3D
,	const float sigma
){
	if(k < 0)
		prms.k = (sigma > 30) ? 8 : 7;
	else
		prms.k = k;

	if(Nf < 0)
		prms.Nf = 4;
	else
		prms.Nf = Nf;

	if(Ns < 0)
		prms.Ns = 7;
	else
		prms.Ns = Ns;
	if(Npr < 0)
		prms.Npr = 5;
	else
		prms.Npr = Npr;
	if(Nb < 0)
		prms.Nb = 2;
	else
		prms.Nb = Nb;

	if(d < 0)
		prms.d = (3*3)/(255.*255.);
	else
		prms.d = (d*d)/(255.*255.);

	if(p < 0)
		prms.p = 4;
	else
		prms.p = p;

	if(N < 0)
		prms.N = 8;
	else
		prms.N = N;

	if(tau < 0)
		prms.tau = (sigma > 30) ? 3000/(255.*255.) : 1500/(255.*255.);
	else
		prms.tau = tau;

	if(T_2D == NONE)
		prms.T_2D = DCT;
	else
		prms.T_2D = T_2D;

	if(T_3D == NONE)
		prms.T_3D = HAAR;
	else
		prms.T_3D = T_3D;
}


/**
 * @file   main.cpp
 * @brief  Main executable file. Do not use lib_fftw to
 *         process DCT.
 *
 * @author MARC LEBRUN  <marc.lebrun@cmla.ens-cachan.fr>
 */


int main(int argc, char **argv)
{
    //! Check if there is the right call for the algorithm

	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.tiff" , "> denoised sequence");

	//! General parameters
	const float fSigma = clo_option("-sigma", 0.f, "Add noise of standard deviation sigma");

	//! VBM3D parameters
	const int kHard = clo_option("-kHard", -1 , "< ");
	const int NfHard = clo_option("-NfHard", -1 , "< ");
	const int NsHard = clo_option("-NsHard", -1 , "< ");
	const int NprHard = clo_option("-NprHard", -1 , "< ");
	const int NbHard = clo_option("-NbHard", -1 , "< ");
	const int pHard = clo_option("-pHard", -1 , "< ");
	const int NHard = clo_option("-NHard", -1 , "< ");
	const int dHard = clo_option("-dHard", -1 , "< ");
	const float tauHard = clo_option("-tauHard", -1. , "< ");
	const int kWien = clo_option("-kWien", -1 , "< ");
	const int NfWien = clo_option("-NfWien", -1 , "< ");
	const int NsWien = clo_option("-NsWien", -1 , "< ");
	const int NprWien = clo_option("-NprWien", -1 , "< ");
	const int NbWien = clo_option("-NbWien", -1 , "< ");
	const int pWien = clo_option("-pWien", -1 , "< ");
	const int NWien = clo_option("-NWien", -1 , "< ");
	const int dWien = clo_option("-dWien", -1 , "< ");
	const float tauWien = clo_option("-tauWien", -1. , "< ");

	const float lambda3D = clo_option("-lambda3d", -1. , "< ");
	const unsigned color_space  =  (unsigned) clo_option("-color", 0 , "< color space");

	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input images.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	//! Variables initialization
	const unsigned T_2D_hard  = (unsigned) clo_option("-T2dh", NONE , "< tau_2D_hard");
	if (T_2D_hard != NONE && T_2D_hard != DCT && T_2D_hard != BIOR)
	{
		cout << "T_2d_hard is not known. Choice is :" << endl;
		cout << " -dct (" << DCT << ")" << endl;
		cout << " -bior (" << BIOR << ")" << endl;
		return EXIT_FAILURE;
	}

	const unsigned T_2D_wien  = (unsigned) clo_option("-T2dw", NONE , "< tau_2D_wien");
	if (T_2D_wien != NONE && T_2D_wien != DCT && T_2D_wien != BIOR)
	{
		cout << "T_2d_wien is not known. Choice is :" << endl;
		cout << " -dct (" << DCT << ")" << endl;
		cout << " -bior (" << BIOR << ")" << endl;
		return EXIT_FAILURE;
	};

	const unsigned T_3D_hard  = (unsigned) clo_option("-T3dh", NONE , "< tau_3D_hard");
	if (T_3D_hard != NONE && T_3D_hard != HAAR && T_3D_hard != HADAMARD)
	{
		cout << "T_3d_hard is not known. Choice is :" << endl;
		cout << " -haar (" << HAAR << ")" << endl;
		cout << " -hadamard (" << HADAMARD << ")" << endl;
		return EXIT_FAILURE;
	}

	const unsigned T_3D_wien  = (unsigned) clo_option("-T3dw", NONE , "< tau_3D_wien");
	if (T_3D_wien != NONE && T_3D_wien != HAAR && T_3D_wien != HADAMARD)
	{
		cout << "T_3d_wien is not known. Choice is :" << endl;
		cout << " -haar (" << HAAR << ")" << endl;
		cout << " -hadamard (" << HADAMARD << ")" << endl;
		return EXIT_FAILURE;
	};

	Parameters prms_1;
	Parameters prms_2;

	initializeParameters_1(prms_1, kHard, NfHard, NsHard, NprHard, NbHard, pHard, NHard, dHard, tauHard, lambda3D, T_2D_hard, T_3D_hard, fSigma);
	initializeParameters_2(prms_2, kWien, NfWien, NsWien, NprWien, NbWien, pWien, NWien, dWien, tauWien, T_2D_wien, T_3D_wien, fSigma);

    if(prms_1.k <= 0)
    {
        cout << "Invalid patch size for the first step" << endl;
        return EXIT_FAILURE;
    }

    // Force both buffer to have the same size for now
    prms_2.Nf = prms_1.Nf;


	//! Preprocessing (KaiserWindow, Threshold, DCT normalization, ...)
    int kHard_2 = prms_1.k*prms_1.k;
	vector<float> kaiser_window_1(kHard_2);
	vector<float> coef_norm_1(kHard_2);
	vector<float> coef_norm_inv_1(kHard_2);
	preProcess(kaiser_window_1, coef_norm_1, coef_norm_inv_1, prms_1.k);

    int kWien_2 = prms_2.k*prms_2.k;
	vector<float> kaiser_window_2(kWien_2);
	vector<float> coef_norm_2(kWien_2);
	vector<float> coef_norm_inv_2(kWien_2);
	preProcess(kaiser_window_2, coef_norm_2, coef_norm_inv_2, prms_2.k);

	//! Preprocessing of Bior table
	vector<float> lpd, hpd, lpr, hpr;
	bior15_coef(lpd, hpd, lpr, hpr);

	//! Declarations
    int w,h,d;
    FILE* in = vpp_init_input(input_path.c_str(), &w, &h, &d);
    if (!in)
        return fprintf(stderr, "cannot initialize input '%s'\n", input_path.c_str()), 1;
    FILE* out = vpp_init_output(final_path.c_str(), w, h, d);
    if (!out)
        return fprintf(stderr, "cannot initialize output '%s'\n", final_path.c_str()), 1;
    
    vector<float*> buffer_input(prms_1.Nf, NULL);
    vector<float*> buffer_basic(prms_2.Nf, NULL);
    float* final_estimate;

    //! Initialize the buffers
    for(int i = 0; i < prms_1.Nf; ++i)
    {
        buffer_input[i] = (float*) malloc(w*h*d*sizeof(*(buffer_input[i])));
        buffer_basic[i] = (float*) malloc(w*h*d*sizeof(*(buffer_basic[i])));
    }
    final_estimate = (float*) malloc(w*h*d*sizeof*final_estimate);

    //! Process the frames
    int size_buffer = 0;
    int index = 0;
    while (vpp_read_frame(in, buffer_input[index], w, h, d)) {
        size_buffer = std::min(size_buffer+1, (int)prms_1.Nf);

        // Change colorspace (RGB to OPP)
        if(color_space == 0)
            transformColorSpace(buffer_input[index], w, h, d, true);

        if (run_vbm3d(fSigma, buffer_input, buffer_basic, final_estimate, w, h, d, prms_1, prms_2, index, size_buffer, kaiser_window_1, coef_norm_1, coef_norm_inv_1, kaiser_window_2, coef_norm_2, coef_norm_inv_2, lpd, hpd, lpr, hpr)
                != EXIT_SUCCESS)
            return EXIT_FAILURE;

        // Inverse the colorspace for the output (OPP to RGB)
        if(color_space == 0)
            transformColorSpace(final_estimate, w, h, d, false);

        // send the frame to the next step
        if (!vpp_write_frame(out, final_estimate, w, h, d))
            break;

        index = (index+1) % prms_1.Nf;
    }

    for(int i = 0; i < prms_1.Nf; ++i)
    {
        free(buffer_input[i]);
        free(buffer_basic[i]);
    }
    free(final_estimate);
    fclose(in);
    fclose(out);

	return EXIT_SUCCESS;
}
