#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "uniform&gaussianRNG.h"
#include "khanh_1112proposal_modulation_code.h"

/*******************  Constant values definition  **************************/
const double PI = 3.14159;                       //// value of π
const double DELTA = 0.000005;                   //// step size that controls the coefficient-update resolution of the equaliser
//const double DELTA          =   0.000005;      //// step size that controls the coefficient-update resolution of the equaliser
												 //// A larger value speeds up adaptation but lowers accuracy; a smaller value gives slower adaptation but higher accuracy.
const int Training_Block = 3;                    //// number of blocks reserved for training
const int M_ary = 4;                             //// total number of symbols (2-ary = 1 bit, 4-ary = 2 bits, 3-ary = 1.5 bits)
const int BLOCK = 5;                           //// total number of blocks => 1000 block
const int Page_Size = 1024;                      //// upper limit of block length, constant
const int Data_Size = Page_Size;                 //// data length per page
const double SigmaB = 1;                         //// standard deviation of block noise => blur noise
const int SELL = 5;                              //// pixel centre = 0; SELL must be odd
const int Unit_Sell = 11;                        //// number of sub-cells in one sell
const int Equalizer_Length = 5;                  //// equaliser length (must be odd)
const int Pridictor_Length = 3;                  //// predictor length
const int Adaptive_Length = 5;                   //// adaptive filter length
const double Mis_const = 1;                      //// constant used to limit mis-alignment
//// 0: random; 1: constant value (Mis_alignment_x, Mis_alignment_y)
const double Mis_alignment_dB = 0; // 0, 10, 20, 30
const double Mis_alignment_x = 0.0;
const double Mis_alignment_y = 0.0;
const int Target_Order = 3;                    //// order of PR target (must be odd)
const int depth = 29;                   //// Viterbi constraint length
const int noise_free = 1;                    //// 0: noise-free; 1: noise present
const double START_SNR = 12;
const double END_SNR = 12;
const double SNR_INTERVAL = 1;
const char* mis[3] = { "mis0", "mis10", "mis20" };
const int checkHist = 0;                         //// 1: gather histogram
const int Encoder = 11;
//// 0: random; 10: 24 mc_chi; 11: 46 mc_chi
const char* s_encoder[10] = { "ran","prof_park","kim_3_4","sun_park_69",
							  "kim_2_3","seungmin_69",
							  "gg_grad_23_modulation_kim",
							  "gg_grad_23_modulation_prof",
							  "gg_dmin2","gg_dmin2_state" };
const int DRAWING_BER_CURVE = 1;                 //// 0: BER in float; 1: BER as integer
double POWER = 1;                                //// normalisation power
const int DETECTOR = 0;

//// 0: 1D detector,

const int Modulation_Decoder = 0;	//// For ECC 4/6 codes. 0: Viterbi, 1: SOVA
const int Iteration_Det_Dem = 1;	//// # of the iteration between detecter(2D SOVA) and demodulation (SOVA for demodulation)
const char* s_detecter[8] = { "1D_Viterbi","2D_NPSOVA","2D NPSOVA Hard","Only AWGN","Adaptive 2D NPSOVA","proposed iterative 2D SOVA","2D Log-MAP","gg_grad_A" };
const int Quantization = 0;								//// 1,2,3,4 : Quantization, else: None Quantization, 

const int NoiseFilter = 0;								//// 0: No noisefilter,  1: use noisefilter
const int Adaptive_yesornot = 0;								//// 0: No adaptive,  1: use adaptive
const int AdaptiveNF_yesornot = 0;								//// 0: No adaptive Nosefilter,  1: use adaptive Nosefilter

const int M_number = 4;								//// 0: no adapt M-alorithm, 1~ : adapt M-alorithm number
const int Constant_period = 0;								//// In M-algirhtm, not use M-algorithm from Constant_period to 0 reversly
const int lee_proposed_method = 0;	//// 0: 1st-8state, 1: 2nd-8state, 2: 3rd-0state, 3: 4th-16state
const int lee_proposed_constraint = 30;	//// 30, 100, 500, 2000

const int SCALE_DEMODUL = 10;								//// During de-modulating, the number applied Viterbi decoder

long k_noise = 1;

//////////////////////// Removed the 03- and 30-state cases from Viterbi //////////////////////////////
const int gg_pr_viterbi = 0;							///     1 : proposed, else : conventional
const int gg_reduce_last = 0;
int ggg = 0;
const int modul_method = 1;							///		0 : state 8, 1 : state 16



void TwoD_Quantization(double** data, double** outdata, int Page_Size, double max, double min, int quantization_number);
void HoloChannel(double** h, double SigmaB, double Mis_alignment_sigma, int SELL, int Unit_Sell);
void PR_Target_set(double** Target_Coef, int target_order);
void PRML_Convolution(double** trellis, double* Target_Coef, int Target_Order);
void Hologram_ChannelOut(double** outCH, double** h, int** inCH, int Page_Size, int SELL, double sigma);
void Hologram_EqualizerOut_x(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order);
void Hologram_EqualizerOut_y(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order);
void NPML_Pridictor_x(double* Pridictor_Coef_x, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length);
void NPML_Pridictor_y(double* Pridictor_Coef_y, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length);
void Hologram_SOVA_x(double** out_SOVA, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth);
void Hologram_SOVA_y(double** out_SOVA, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth);
void MakeInterleaverFile(int* Int, int* DeInt, int nSize);
void Interleaver_for_Holo(int** Page, int** Inter_Page, int* Inter_table, int Page_Size);
void Interleaver_for_Holo(double** Page, double** Inter_Page, int* Inter_table, int Page_Size);

void PRML_Convolution_gg_grad_23(double** trellis, double* Target_Coef, int Target_Order, int input_size, int* input);
void PRML_Convolution_for_modul_gg_grad_23(int*** trellis_modul);
void PRML_Convolution_for_modul_gg_grad_23_prof(int*** trellis_modul);
void PRML_Convolution_for_modul_gg_grad_23_prof_method0(int**** trellis_modul);
void PRML_Convolution_for_modul_gg_dmin2(int*** trellis_modul);
void PRML_Convolution_for_modul_gg_dmin2_method0(int**** trellis_modul);
void Hologram_SOVA_x_for_gg_grad(double** out_SOVA, double** eq_out, double** EI, double** trellis, double** trellisA, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth);



void main()
{
	clock_t c_time;
	c_time = clock(); // c_time Check Start 

	int i, j, k;

	int A_input[2] = { 0,1 };
	int BC_input[4] = { 0,1,2,3 };
	int total_states_method0 = 0;

	double** trellis_A;


	int file_see = 0;

	FILE* OUT, * BRIEF_RESULT, * BRIEF_RESULT_MIS, * fber, * CHOUT[M_ary], * gg_fp_en, * gg_fp_de, * gg_fp_pattern;

	fber = fopen("test.dat", "w");
	BRIEF_RESULT = fopen("breif_result.m", "a");
	BRIEF_RESULT_MIS = fopen("breif_result_mis.m", "a");
	gg_fp_en = fopen("encoded.dat", "a");
	gg_fp_de = fopen("decoded.dat", "a");
	gg_fp_pattern = fopen("pattern.dat", "a");

	char nameMis0[255];
	char nameMis1[255];
	char nameMis2[255];
	char nameMis3[255];
	char nameMis4[255];
	char nameMis5[255];

	FILE* fber_outCH, * fber_outCH_Mary;
	FILE* fber_outSOVA, * fber_outSOVA_Mary;
	FILE* fber_outDEM, * fber_outDEM_Mary;

	if (Mis_alignment_x == 0)
	{
		sprintf(nameMis0, "%s_fber_outCH.txt", mis[0]);
		sprintf(nameMis1, "%s_fber_outCH_Mary.txt", mis[0]);
		sprintf(nameMis2, "%s_fber_outSOVA.txt", mis[0]);
		sprintf(nameMis3, "%s_fber_outSOVA_Mary.txt", mis[0]);
		sprintf(nameMis4, "%s_fber_outDEM.txt", mis[0]);
		sprintf(nameMis5, "%s_fber_outDEM_Mary.txt", mis[0]);

		fber_outCH = fopen(nameMis0, "a");
		fber_outCH_Mary = fopen(nameMis1, "a");

		fber_outSOVA = fopen(nameMis2, "a");
		fber_outSOVA_Mary = fopen(nameMis3, "a");

		fber_outDEM = fopen(nameMis4, "a");
		fber_outDEM_Mary = fopen(nameMis5, "a");

	}
	else if (Mis_alignment_x == 0.1)
	{
		sprintf(nameMis0, "%s_fber_outCH.txt", mis[1]);
		sprintf(nameMis1, "%s_fber_outCH_Mary.txt", mis[1]);
		sprintf(nameMis2, "%s_fber_outSOVA.txt", mis[1]);
		sprintf(nameMis3, "%s_fber_outSOVA_Mary.txt", mis[1]);
		sprintf(nameMis4, "%s_fber_outDEM.txt", mis[1]);
		sprintf(nameMis5, "%s_fber_outDEM_Mary.txt", mis[1]);

		fber_outCH = fopen(nameMis0, "a");
		fber_outCH_Mary = fopen(nameMis1, "a");

		fber_outSOVA = fopen(nameMis2, "a");
		fber_outSOVA_Mary = fopen(nameMis3, "a");

		fber_outDEM = fopen(nameMis4, "a");
		fber_outDEM_Mary = fopen(nameMis5, "a");
	}
	else
	{
		sprintf(nameMis0, "%s_fber_outCH.txt", mis[2]);
		sprintf(nameMis1, "%s_fber_outCH_Mary.txt", mis[2]);
		sprintf(nameMis2, "%s_fber_outSOVA.txt", mis[2]);
		sprintf(nameMis3, "%s_fber_outSOVA_Mary.txt", mis[2]);
		sprintf(nameMis4, "%s_fber_outDEM.txt", mis[2]);
		sprintf(nameMis5, "%s_fber_outDEM_Mary.txt", mis[2]);

		fber_outCH = fopen(nameMis0, "a");
		fber_outCH_Mary = fopen(nameMis1, "a");

		fber_outSOVA = fopen(nameMis2, "a");
		fber_outSOVA_Mary = fopen(nameMis3, "a");

		fber_outDEM = fopen(nameMis4, "a");
		fber_outDEM_Mary = fopen(nameMis5, "a");
	}






	unsigned long count = 0;
	unsigned long block_count = 0;
	unsigned long count_M_ary = 0;
	unsigned long count_modul = 0;
	unsigned long count_modul_M_ary = 0;
	unsigned long error = 0;
	unsigned long error_M_ary = 0;
	unsigned long error_eq = 0;
	unsigned long error_eq_other = 0;
	unsigned long error_viterbi[Iteration_Det_Dem * 2] = { 0, };
	unsigned long error_viterbi_M_ary[Iteration_Det_Dem * 2] = { 0, };
	unsigned long error_viterbi_M_ary_gg[6] = { 0, };
	unsigned long error_demodul[Iteration_Det_Dem * 2] = { 0, };
	unsigned long error_demodul_M_ary[Iteration_Det_Dem * 2] = { 0, };
	unsigned long error_demodul_check[Iteration_Det_Dem] = { 0, };

	int temp_count[M_ary] = { 0, };
	double temp_sum[M_ary] = { 0, };

	double mis_sigma = 0.0;
	double sigma = 0.0;										//// SNR_DB¸¦ À̀¿ëÇØ ±¸ÇÑ sigma¸¦ ÀúÀåÇÏ´Â °÷
	double snr;												//// SNR_DB¸¦ À̀¿ëÇØ SNRÀÇ °ªÀ» ÀúÀåÇÏ´Â °÷
	double eqout_total_x = 0, eqout_average_x = 0, eqout_total_y = 0, eqout_average_y = 0;
	double chout_average[M_ary - 1] = { 0, };
	char strH[255];

	double** h;												//// Ă¤³ÎÀ» ÀúÀå
	double** Target_Coef;									//// PR TargetÀÇ coefficient¸¦ Á¤ÇÑ´Ù.
	double** trellis;

	Target_Coef = (double**)malloc(Target_Order * sizeof(double*));
	for (i = 0; i < Target_Order; i++)
		Target_Coef[i] = (double*)malloc(Target_Order * sizeof(double));

	////// For Viterbi and BCJR ///////
	trellis = (double**)malloc((int)pow((double)M_ary, Target_Order - 1) * sizeof(double*));
	for (i = 0; i < pow((double)M_ary, Target_Order - 1); i++)
		trellis[i] = (double*)calloc(M_ary * 2, sizeof(double));

	////// For Viterbi and BCJR : gukhui's graduation ///////
	trellis_A = (double**)malloc((int)pow((double)M_ary - 2, Target_Order - 1) * sizeof(double*));
	for (i = 0; i < pow((double)M_ary - 2, Target_Order - 1); i++)
		trellis_A[i] = (double*)calloc((M_ary - 2) * 2, sizeof(double));

	int*** trellis_modul = (int***)malloc((int)pow((double)M_ary, 2) * sizeof(int**));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		trellis_modul[i] = (int**)malloc((int)pow((double)M_ary, 2) * sizeof(int*));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			trellis_modul[i][j] = (int*)calloc(3, sizeof(int));

	int*** trellis_modul_prof = (int***)malloc((int)pow((double)M_ary, 2) * sizeof(int**));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		trellis_modul_prof[i] = (int**)malloc((int)pow((double)M_ary, 2) * sizeof(int*));
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			trellis_modul_prof[i][j] = (int*)calloc(3, sizeof(int));

	total_states_method0 = ((int)pow((double)M_ary, 2)) / 2;
	int**** trellis_modul_prof_method0 = (int****)malloc(total_states_method0 * sizeof(int***));
	for (i = 0; i < total_states_method0; i++)
		trellis_modul_prof_method0[i] = (int***)malloc(total_states_method0 * sizeof(int**));
	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			trellis_modul_prof_method0[i][j] = (int**)malloc(2 * sizeof(int*));
	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			for (k = 0; k < 2; k++)
				trellis_modul_prof_method0[i][j][k] = (int*)calloc(3, sizeof(int));


	///// For channel modeling ///////
	h = (double**)malloc(SELL * sizeof(double*));
	for (i = 0; i < SELL; i++)
		h[i] = (double*)calloc(SELL, sizeof(double));

	///// For PR equalizer //////
	double** Equalizer_Coef_x = (double**)malloc(Equalizer_Length * sizeof(double*));	//// EqualizerÀÇ coefficient¸¦ ÀúÀåÇÏ´Â °÷
	double** Equalizer_Coef_y = (double**)malloc(Equalizer_Length * sizeof(double*));	//// EqualizerÀÇ coefficient¸¦ ÀúÀåÇÏ´Â °÷
	for (i = 0; i < Equalizer_Length; i++)
		Equalizer_Coef_x[i] = (double*)calloc(Equalizer_Length, sizeof(double)), Equalizer_Coef_y[i] = (double*)calloc(Equalizer_Length, sizeof(double));

	double* Pridictor_Coef_x = (double*)calloc(Pridictor_Length + 1, sizeof(double));	//// ÇÁ¸®µñÅÍÀÇ coefficient¸¦ ÀúÀåÇÏ´Â °÷
	double* Pridictor_Coef_y = (double*)calloc(Pridictor_Length + 1, sizeof(double));	//// ÇÁ¸®µñÅÍÀÇ coefficient¸¦ ÀúÀåÇÏ´Â °÷

	///// For modulation codes /////
	int* input_park_2_3_4ary = (int*)calloc(Data_Size * (Data_Size / 3) * 2, sizeof(int));
	int* output_park_2_3_4ary = (int*)calloc(Data_Size * (Data_Size / 3) * 2, sizeof(int));
	int* input_kim_3_4_4ary = (int*)calloc((Data_Size / 2) * (Data_Size / 2) * 3, sizeof(int));
	int* output_kim_3_4_4ary = (int*)calloc((Data_Size / 2) * (Data_Size / 2) * 3, sizeof(int));
	int* input_sun_park_6_9_4ary = (int*)calloc((Data_Size / 3) * (Data_Size / 3) * 6, sizeof(int));
	int* output_sun_park_6_9_4ary = (int*)calloc((Data_Size / 3) * (Data_Size / 3) * 6, sizeof(int));

	int* input_khanh_1112_4ary = (int*)calloc((Data_Size / 4) * (Data_Size / 3) * 11, sizeof(int));
	int* output_khanh_1112_4ary = (int*)calloc((Data_Size / 4) * (Data_Size / 3) * 11, sizeof(int));

	int* input_chi_24_4ary = (int*)calloc(Data_Size * (Data_Size / 4) * 2, sizeof(int));
	int* output_chi_24_4ary = (int*)calloc(Data_Size * (Data_Size / 4) * 2, sizeof(int));


	double** EI = (double**)malloc(Page_Size * sizeof(double*));
	int** outLDPC = (int**)malloc(Page_Size * sizeof(int*));
	for (i = 0; i < Page_Size; i++)
		outLDPC[i] = (int*)calloc(Page_Size, sizeof(int)), EI[i] = (double*)calloc(Page_Size, sizeof(double));

	int** inLDPC = (int**)malloc(Data_Size * sizeof(int*));							//// binary input
	for (i = 0; i < Data_Size; i++)
		inLDPC[i] = (int*)calloc(Data_Size * M_ary, sizeof(int));

	int** inLDPC_symbol = (int**)malloc(Data_Size * sizeof(int*));							//// M-ary symbol input
	for (i = 0; i < Data_Size; i++)
		inLDPC_symbol[i] = (int*)calloc(Data_Size, sizeof(int));

	int** inCH = (int**)malloc((Page_Size + (Equalizer_Length / 2) * 2) * sizeof(int*));	//// Ă¤³ÎÀ» Åë°úÇ̉ µ¥À̀ÅÍ
	for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
		inCH[i] = (int*)calloc(Page_Size + (Equalizer_Length / 2) * 2, sizeof(int));

	double** outCH = (double**)malloc((Page_Size + (Equalizer_Length / 2) * 2) * sizeof(double*));											//// Ă¤³ÎÀ» Åë°úÇÑ ÈÄÀÇ µ¥À̀ÅÍ
	for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
		outCH[i] = (double*)calloc((Page_Size + (Equalizer_Length / 2) * 2), sizeof(double));

	double** outEQ_x = (double**)malloc(Page_Size * sizeof(double*));
	double** outEQ_y = (double**)malloc(Page_Size * sizeof(double*));
	int** n_outViterbi = (int**)malloc(Page_Size * sizeof(int*));
	double** d_outViterbi = (double**)malloc(Page_Size * sizeof(double));
	double** d_outViterbi_x = (double**)malloc(Page_Size * sizeof(double));
	double** d_outViterbi_y = (double**)malloc(Page_Size * sizeof(double));

	for (i = 0; i < Page_Size; i++)
	{
		outEQ_x[i] = (double*)calloc(Page_Size, sizeof(double));
		outEQ_y[i] = (double*)calloc(Page_Size, sizeof(double));
		n_outViterbi[i] = (int*)calloc(Page_Size, sizeof(int));
		d_outViterbi[i] = (double*)calloc(Page_Size, sizeof(double));
		d_outViterbi_x[i] = (double*)calloc(Page_Size, sizeof(double));
		d_outViterbi_y[i] = (double*)calloc(Page_Size, sizeof(double));

	}

	int* Interleaver_table = (int*)malloc(Page_Size * Page_Size * sizeof(int));
	int* De_Interleaver_table = (int*)malloc(Page_Size * Page_Size * sizeof(int));

	double** d_Interleaved = (double**)malloc(Page_Size * sizeof(double*));
	int** n_Interleaved = (int**)malloc(Page_Size * sizeof(int*));
	for (i = 0; i < Page_Size; i++)
	{
		d_Interleaved[i] = (double*)calloc(Page_Size, sizeof(double));
		n_Interleaved[i] = (int*)calloc(Page_Size, sizeof(int));
	}

	int** en_piece_sun_park;
	int** de_piece_sun_park;
	if (Encoder == 3)
	{
		en_piece_sun_park = (int**)malloc(64 * sizeof(int*));
		de_piece_sun_park = (int**)malloc(64 * sizeof(int*));
		for (i = 0; i < 64; i++)
		{
			en_piece_sun_park[i] = (int*)calloc(4, sizeof(int));
			de_piece_sun_park[i] = (int*)calloc(3, sizeof(int));
		}
		for (i = 0; i < 64; i++)
		{
			for (j = 0; j < 4; j++)
				en_piece_sun_park[i][j] = (i / int(pow(3.0, 3 - j))) % 3;
			for (j = 0; j < 3; j++)
				de_piece_sun_park[i][j] = (i >> ((2 - j) * 2)) % 4;
		}
	}

	FILE* OUTchi_24;


	int** piece_chi_24 = (int**)malloc(16 * sizeof(int*));
	for (i = 0; i < 16; i++)
		piece_chi_24[i] = (int*)calloc(4, sizeof(int));
	if (10 == Encoder)
	{
		OUTchi_24 = fopen("4ary24mc_16.txt", "r");
		for (i = 0; i < 16; i++)
			for (j = 0; j < 4; j++)
				fscanf(OUTchi_24, "%d", &piece_chi_24[i][j]);
		fclose(OUTchi_24);
	}
	//MakeInterleaverFile(Interleaver_table, De_Interleaver_table, Page_Size/3*Page_Size/2);

	snr = pow(10, (Mis_alignment_dB / 10.0));
	mis_sigma = sqrt(float(POWER / (snr)));

	HoloChannel(h, SigmaB, mis_sigma, SELL, Unit_Sell);
	PR_Target_set(Target_Coef, Target_Order);
	PRML_Convolution(trellis, Target_Coef[Target_Order / 2], Target_Order);

	if (Encoder == 6)
		PRML_Convolution_for_modul_gg_grad_23(trellis_modul);
	if (Encoder == 7)
	{
		PRML_Convolution_for_modul_gg_grad_23_prof(trellis_modul_prof);
		PRML_Convolution_for_modul_gg_grad_23_prof_method0(trellis_modul_prof_method0);
	}
	if (Encoder == 9)
	{
		PRML_Convolution_for_modul_gg_dmin2(trellis_modul_prof);
		PRML_Convolution_for_modul_gg_dmin2_method0(trellis_modul_prof_method0);
	}

	PRML_Convolution_gg_grad_23(trellis_A, Target_Coef[Target_Order / 2], Target_Order, 2, A_input);

	printf("Multilevel Hologram sigmaB %.2lf, Sell %d, unit_sell %d, Target_Order %d, Equalizer_Length %d\n", SigmaB, SELL, Unit_Sell, Target_Order, Equalizer_Length);
	printf("Quantization %d, page_size %d, block %d, total_bit %d\n", Quantization, Page_Size, BLOCK, long(Page_Size * Page_Size * BLOCK));
	printf("Adaptive_Length %d, Pridictor_Length %d, NoiseFilter %d, Adaptive %d, NF Adaptive %d\nDELTA %lf, ", Adaptive_Length, Pridictor_Length, NoiseFilter, Adaptive_yesornot, AdaptiveNF_yesornot, DELTA);
	printf("Detector is %s, ", s_detecter[DETECTOR]);
	if (Mis_const == 0)
		printf("Random Misalignment %lf db, \n", Mis_alignment_dB);
	else if (Mis_const == 1)
		printf("Misalignment (%1.2lf %1.2lf), \n", Mis_alignment_x, Mis_alignment_y);
	printf("Target is\n");
	for (i = 0; i < Target_Order; i++)
	{
		for (j = 0; j < Target_Order; j++)
			printf("%0.1lf ", Target_Coef[i][j]);
		printf("\n");
	}
	fprintf(fber, "Multilevel Hologram sigmaB %.2lf, Sell %d, unit_sell %d, Target_Order %d, Equalizer_Length %d\n", SigmaB, SELL, Unit_Sell, Target_Order, Equalizer_Length);
	fprintf(fber, "Quantization %d, page_size %d, block %d, total_bit %d\n", Quantization, Page_Size, BLOCK, long(Page_Size * Page_Size * BLOCK));
	fprintf(fber, "Adaptive_Length %d, Pridictor_Length %d, NoiseFilter %d, Adaptive %d, NF Adaptive %d\nDELTA %lf, ", Adaptive_Length, Pridictor_Length, NoiseFilter, Adaptive_yesornot, AdaptiveNF_yesornot, DELTA);
	fprintf(fber, "Detector is %s, ", s_detecter[DETECTOR]);
	if (Mis_const == 0)
		fprintf(fber, "Random Misalignment %lf db, \n", Mis_alignment_dB);
	else if (Mis_const == 1)
		fprintf(fber, "Misalignment (%1.2lf %1.2lf), \n", Mis_alignment_x, Mis_alignment_y);
	fprintf(fber, "Target is\n");
	for (i = 0; i < Target_Order; i++)
	{
		for (j = 0; j < Target_Order; j++)
			fprintf(fber, "%0.1lf ", Target_Coef[i][j]);
		fprintf(fber, "\n");
	}

	if (Encoder == 0)
	{
		printf("No modulation (random) code\n");
		fprintf(fber, "No modulation (random) code\n");
	}
	else if (Encoder == 1)
	{
		printf("Prof. Park 2/3 modulation code\n");
		fprintf(fber, "Prof. Park 2/3 modulation code\n");
	}
	else if (Encoder == 2)
	{
		printf("Kim's simple 3/4 modulation code\n");
		fprintf(fber, "Kim's simple 3/4 modulation code\n");
	}
	else if (Encoder == 3)
	{
		printf("Sunsuk & Gunhwan's 6/9 modulation code\n");
		fprintf(fber, "Sunsuk & Gunhwan's 6/9 modulation code\n");
	}
	else if (Encoder == 4)
	{
		printf("Kim odd even 2/3 modulation code\n");
		fprintf(fber, "Kim odd even 2/3 modulation code\n");
	}
	else if (Encoder == 5)
	{
		printf("Seungmin's 6/9 modulation code\n");
		fprintf(fber, "Seungmin's 6/9 modulation code\n");
	}
	else if (Encoder == 11)
	{
		printf("khanh's 11/12 modulation code\n");
		//fprintf(fber, "Seungmin's 6/9 modulation code\n");
	}
	else if (Encoder == 6)
	{
		printf("gg_grad 2/3 modulation code\n");
		fprintf(fber, "gg_grad 2/3 modulation code\n");
	}
	else if (Encoder == 7)
	{
		printf("gg_grad 2/3 modulation code based on prof.park code\n");
		fprintf(fber, "gg_grad 2/3 modulation code based on prof.park code\n");
	}
	else if (Encoder == 8)
	{
		printf("gg_dmin2_code\n");
		fprintf(fber, "gg_dmin2_code\n");
	}
	else if (Encoder == 9)
	{
		printf("gg_dmin2_state_code based on modul method %d\n", modul_method);
		fprintf(fber, "gg_dmin2_state_code based on modul method %d\n", modul_method);
	}
	else if (Encoder == 10)
	{
		printf("Chi 24 proposal\n");
	}
	for (double SNR_DB = START_SNR - SNR_INTERVAL; SNR_DB < END_SNR;)
	{
		SNR_DB += SNR_INTERVAL;

		// Initialize variable
		count = 0, count_M_ary = 0, count_modul = 0, count_modul_M_ary = 0;
		block_count = 0, error = 0, error_eq = 0, error_eq_other = 0, error_M_ary = 0;
		error_viterbi[0] = 0, error_viterbi_M_ary[0] = 0;
		error_viterbi[1] = 0, error_viterbi_M_ary[1] = 0;
		error_demodul[0] = 0, error_demodul_M_ary[0] = 0;
		error_demodul[1] = 0, error_demodul_M_ary[1] = 0;

		for (i = 0; i < Equalizer_Length; i++)
			for (j = 0; j < Equalizer_Length; j++)
				Equalizer_Coef_x[i][j] = 0, Equalizer_Coef_y[i][j] = 0;

		snr = pow(10, (SNR_DB / 10.0));
		sigma = float(POWER / (snr));

		while (block_count < Training_Block && DETECTOR != 3)
		{
			srand(0);
			if (Encoder == 0)
			{
				if (M_ary == 2)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC[i][j] = inLDPC_symbol[i][j] = rand() % 2, count++, count_M_ary++;
				}
				else if (M_ary == 4)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size * 2; j++)
							inLDPC[i][j] = rand() % 2, count++;
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC_symbol[i][j] = (inLDPC[i][j * 2] << 1) + inLDPC[i][j * 2 + 1], count_M_ary++;
				}
			}



			else if (Encoder == 11)
			{
				if (M_ary != 4)
				{
					printf("\n\n mismatch encoder and M_ary\n You should select Encoder 3 and 4-ary.\n\n");
					exit(0);
				}
				for (i = 0; i < (Data_Size / 4) * (Data_Size / 3) * 11; i++)
					input_khanh_1112_4ary[i] = rand() % M_ary, count_modul += 2, count_modul_M_ary++;
				Encode_khanhProposal1112_4ary(inLDPC_symbol, input_khanh_1112_4ary, Page_Size);
				count_M_ary += Page_Size * Page_Size;
				count += Page_Size * Page_Size * 2;
			}

			for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
				for (j = 0; j < Page_Size + (Equalizer_Length / 2) * 2; j++)
					if ((i > Equalizer_Length / 2 - 1 && i < Page_Size + Equalizer_Length / 2) && (j > Equalizer_Length / 2 - 1 && j < Page_Size + Equalizer_Length / 2))
						inCH[i][j] = inLDPC_symbol[i - Equalizer_Length / 2][j - Equalizer_Length / 2];
					else
						inCH[i][j] = 0;

			Hologram_ChannelOut(outCH, h, inCH, Page_Size + (Equalizer_Length / 2) * 2, SELL, sigma);
			// Hologram_EqualizerOut_x(outEQ_x, Equalizer_Coef_x, outCH, Page_Size, SELL, Equalizer_Length, 0, inCH, Target_Coef, Target_Order);
			Hologram_EqualizerOut_y(outEQ_y, Equalizer_Coef_y, outCH, Page_Size, SELL, Equalizer_Length, 0, inCH, Target_Coef, Target_Order);
			if (NoiseFilter == 1)
			{
				NPML_Pridictor_x(Pridictor_Coef_x, outEQ_x, Page_Size, SELL, Equalizer_Length, inCH, Target_Coef, Target_Order, Pridictor_Length);
				NPML_Pridictor_y(Pridictor_Coef_y, outEQ_y, Page_Size, SELL, Equalizer_Length, inCH, Target_Coef, Target_Order, Pridictor_Length);
			}
			block_count++;
		}
		block_count = 0, count = 0, count_M_ary = 0;
		count_modul = 0, count_modul_M_ary = 0;
		while (block_count < BLOCK)
		{
			if ((block_count + 1) % 50 == 0)
				printf(">");
			if ((block_count + 1) % 200 == 0)
				printf("%d", (block_count + 1) / 100);
			srand(block_count);
			if (Encoder == 0)
			{
				if (M_ary == 2)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC[i][j] = inLDPC_symbol[i][j] = rand() % 2, count++, count_M_ary++;
				}
				else if (M_ary == 4)
				{
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size * 2; j++)
							inLDPC[i][j] = rand() % 2, count++;
					for (i = 0; i < Data_Size; i++)
						for (j = 0; j < Data_Size; j++)
							inLDPC_symbol[i][j] = (inLDPC[i][j * 2] << 1) + inLDPC[i][j * 2 + 1], count_M_ary++;
				}
			}


			else if (Encoder == 11)
			{
				if (M_ary != 4)
				{
					printf("\n\n mismatch encoder and M_ary\n You should select Encoder 3 and 4-ary.\n\n");
					exit(0);
				}
				for (i = 0; i < (Data_Size / 4) * (Data_Size / 3) * 11; i++)
					input_khanh_1112_4ary[i] = rand() % M_ary, count_modul += 2, count_modul_M_ary++;
				Encode_khanhProposal1112_4ary(inLDPC_symbol, input_khanh_1112_4ary, Page_Size);
				count_M_ary += Page_Size * Page_Size;
				count += Page_Size * Page_Size * 2;
			}

			if (file_see == 1)
			{
				sprintf(strH, "inLDPC_symbol_%1d.dat", block_count);
				OUT = fopen(strH, "w");
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
						fprintf(OUT, "%d ", inLDPC_symbol[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);
			}

			for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
				for (j = 0; j < Page_Size + (Equalizer_Length / 2) * 2; j++)
					if ((i > Equalizer_Length / 2 - 1 && i < Page_Size + Equalizer_Length / 2) && (j > Equalizer_Length / 2 - 1 && j < Page_Size + Equalizer_Length / 2))
						inCH[i][j] = inLDPC_symbol[i - Equalizer_Length / 2][j - Equalizer_Length / 2];
					else
						inCH[i][j] = 0;

			if (file_see == 1)
			{
				OUT = fopen("inCH.dat", "w");
				for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
				{
					for (j = 0; j < Page_Size + (Equalizer_Length / 2) * 2; j++)
						fprintf(OUT, "%d ", inCH[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);
			}

			/////////////////////////////// channel out ///////////////////////////////////////////////
			HoloChannel(h, SigmaB, mis_sigma, SELL, Unit_Sell);
			Hologram_ChannelOut(outCH, h, inCH, Page_Size + (Equalizer_Length / 2) * 2, SELL, sigma);

			for (i = 0; i < M_ary; i++)
				temp_count[i] = 0, temp_sum[i] = 0;
			for (i = Equalizer_Length / 2; i < Page_Size + (Equalizer_Length / 2); i++)
			{
				for (j = Equalizer_Length / 2; j < Page_Size + (Equalizer_Length / 2); j++)
					for (k = 0; k < M_ary; k++)
						if (inLDPC_symbol[i - Equalizer_Length / 2][j - Equalizer_Length / 2] == k)
							temp_sum[k] += outCH[i][j], temp_count[k]++;
			}
			for (k = 0; k < M_ary - 1; k++)
				chout_average[k] = ((temp_sum[k] / temp_count[k]) + (temp_sum[k + 1] / temp_count[k + 1])) / 2;
			//chout_average=(temp_sum_0+temp_sum_1)/(temp_count_0+temp_count_1);

			if (checkHist == 1)
			{
				for (k = 0; k < M_ary; k++)
				{
					sprintf(strH, "chout_%d.dat", k);
					CHOUT[k] = fopen(strH, "a");
				}
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
					{
						for (k = 0; k < M_ary; k++)
						{
							if (inLDPC_symbol[i][j] == k)
							{
								fprintf(CHOUT[k], "%.10lf \n", outCH[i + Equalizer_Length / 2][j + Equalizer_Length / 2] / 10);
								break;
							}

						}
					}
				}
				for (k = 0; k < M_ary; k++)
					fclose(CHOUT[k]);
			}

			/////////////////////////////// equalization ///////////////////////////////////////////////
		//	Hologram_EqualizerOut_x(outEQ_x, Equalizer_Coef_x, outCH, Page_Size, SELL, Equalizer_Length, 1, inCH, Target_Coef, Target_Order);
			Hologram_EqualizerOut_y(outEQ_y, Equalizer_Coef_y, outCH, Page_Size, SELL, Equalizer_Length, 1, inCH, Target_Coef, Target_Order);

			if (file_see == 1)
			{
				for (k = 0; k < M_ary; k++)
				{
					sprintf(strH, "eqout_%d.dat", k);
					CHOUT[k] = fopen(strH, "a");
				}
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
					{
						for (k = 0; k < M_ary; k++)
						{
							if (inLDPC_symbol[i][j] == k)
							{
								fprintf(CHOUT[k], "%.10lf \n", outEQ_y[i][j]);
								break;
							}

						}
					}
				}
				for (k = 0; k < M_ary; k++)
					fclose(CHOUT[k]);
			}

			/////////////////////////////// detection ///////////////////////////////////////////////
			if (DETECTOR == 0)
			{
		//		Hologram_SOVA_x(d_outViterbi_x, outEQ_x, EI, trellis, Pridictor_Coef_x, Pridictor_Length, Target_Order, Page_Size, depth);
				Hologram_SOVA_y(d_outViterbi_y, outEQ_y, EI, trellis, Pridictor_Coef_y, Pridictor_Length, Target_Order, Page_Size, depth);

			}
			else if (DETECTOR == 7)
			{
				Hologram_SOVA_x_for_gg_grad(d_outViterbi_x, outEQ_x, EI, trellis, trellis_A, Pridictor_Coef_x, Pridictor_Length, Target_Order, Page_Size, depth);
				ggg = 0;
			}

			if (file_see == 1)
			{
				OUT = fopen("d_outViterbi_x.dat", "w");
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
						fprintf(OUT, "%lf ", d_outViterbi_x[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);

				OUT = fopen("d_outViterbi_y.dat", "w");
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < Page_Size; j++)
						fprintf(OUT, "%lf ", d_outViterbi_y[i][j]);
					fprintf(OUT, "\n");
				}
				fclose(OUT);
			}

			/////////////////////////////// demodulation //////////////////////////////////////////////	

			if (Encoder == 11)
			{
				/*Decode_chiProposal46_4ary(output_chi_46_4ary, d_outViterbi_x, Page_Size);
				for (i = 0; i < (Data_Size / 3)*(Data_Size / 2) * 4; i++)
				{
					if (input_chi_46_4ary[i] != output_chi_46_4ary[i])
						error_demodul_M_ary[0]++;
					if (input_chi_46_4ary[i] >> 1 != output_chi_46_4ary[i] >> 1)
						error_demodul[0]++;
					if (input_chi_46_4ary[i] % 2 != output_chi_46_4ary[i] % 2)
						error_demodul[0]++;
				}	*/

				/*Decode_chiProposal46_4ary(output_chi_46_4ary, d_outViterbi_x, Page_Size);
				for (i = 0; i < (Data_Size / 4) * (Data_Size / 3) * 11; i++)
				{
					if (input_chi_46_4ary[i] != output_chi_46_4ary[i])
						error_demodul_M_ary[0]++;
					if (input_chi_46_4ary[i] >> 1 != output_chi_46_4ary[i] >> 1)
						error_demodul[0]++;
					if (input_chi_46_4ary[i] % 2 != output_chi_46_4ary[i] % 2)
						error_demodul[0]++;
				}*/

				Decode_khanhProposal1112_4ary(output_khanh_1112_4ary, d_outViterbi_y, Page_Size);
				for (i = 0; i < (Data_Size / 4) * (Data_Size / 3) * 11; i++)
				{
					if (input_khanh_1112_4ary[i] != output_khanh_1112_4ary[i])
						error_demodul_M_ary[1]++;
					if (input_khanh_1112_4ary[i] >> 1 != output_khanh_1112_4ary[i] >> 1)
						error_demodul[1]++;
					if (input_khanh_1112_4ary[i] % 2 != output_khanh_1112_4ary[i] % 2)
						error_demodul[1]++;
				}

			}


			/////////////////////////////// error check ///////////////////////////////////////////////
			if (M_ary == 4)
			{
				for (i = 0; i < Page_Size; i++)
				{
					for (j = 0; j < (Page_Size); j++)
					{
						ggg = 0;
						for (k = 0; k < M_ary - 1; k++)
						{
							if (outCH[i + Equalizer_Length / 2][j + Equalizer_Length / 2] < chout_average[k])
							{
								temp_count[0] = k;
								break;
							}
							else
								temp_count[0] = M_ary - 1;
						}
						if (inLDPC_symbol[i][j] != temp_count[0])
							error_M_ary++;
						if (inLDPC_symbol[i][j] >> 1 != temp_count[0] >> 1)
							error++;
						if (inLDPC_symbol[i][j] % 2 != temp_count[0] % 2)
							error++;

						/*if (inLDPC_symbol[i][j] != (int)d_outViterbi_x[i][j])
							error_viterbi_M_ary[0]++;
						if (inLDPC_symbol[i][j] >> 1 != (int)d_outViterbi_x[i][j] >> 1)
							error_viterbi[0]++;
						if (inLDPC_symbol[i][j] % 2 != (int)d_outViterbi_x[i][j] % 2)
							error_viterbi[0]++;*/

						if (inLDPC_symbol[i][j] != (int)d_outViterbi_y[i][j])
							error_viterbi_M_ary[1]++;
						if (inLDPC_symbol[i][j] >> 1 != (int)d_outViterbi_y[i][j] >> 1)
							error_viterbi[1]++;
						if (inLDPC_symbol[i][j] % 2 != (int)d_outViterbi_y[i][j] % 2)
							error_viterbi[1]++;

					}

				}
			}

			block_count++;
		}
		printf("\n");
		printf("count : %d\n count_M_ary : %d\n count_modul : %d\n count_modul_M_ary : %d\n\n",count, count_M_ary, count_modul, count_modul_M_ary);
		printf("snr_db : %.1lf\n", SNR_DB);
		printf("snr_db : %.1lf\t\tMis snr_db : %.1lf\n", SNR_DB, Mis_alignment_dB);
		printf("channel error : %d\tBER : %.12lf\n", error, (double)error / count);
		printf("channel M_ary error : %d\tBER : %.12lf\n", error_M_ary, (double)error_M_ary / count_M_ary);
		printf("viterbi error x : %d\tBER : %.12lf\n", error_viterbi[0], (double)error_viterbi[0] / count);
		printf("viterbi M_ary error x : %d\tBER : %.12lf\n", error_viterbi_M_ary[0], (double)error_viterbi_M_ary[0] / count_M_ary);
		printf("viterbi error y : %d\tBER : %.12lf\n", error_viterbi[1], (double)error_viterbi[1] / count);
		printf("viterbi M_ary error y : %d\tBER : %.12lf\n", error_viterbi_M_ary[1], (double)error_viterbi_M_ary[1] / count_M_ary);
		if (Encoder != 0)
		{
			printf("modulation error x : %d\tBER : %.12lf\n",error_demodul[0],(double)error_demodul[0]/count_modul);
			printf("modulation M_ary error x : %d\tBER : %.12lf\n",error_demodul_M_ary[0],(double)error_demodul_M_ary[0]/count_modul_M_ary);
			printf("modulation error y : %d\tBER : %.12lf\n", error_demodul[1], (double)error_demodul[1] / count_modul);
			printf("modulation M_ary error y : %d\tBER : %.12lf\n", error_demodul_M_ary[1], (double)error_demodul_M_ary[1] / count_modul_M_ary);
		}

		fprintf(fber, "\n");
		fprintf(fber, "snr_db : %.1lf\t\tMis snr_db : %.1lf\n", SNR_DB, Mis_alignment_dB);
		fprintf(fber_outCH, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error / count);
		fprintf(fber_outCH_Mary, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_M_ary / count_M_ary);
		fprintf(fber,"viterbi error x : %d\tBER : %.12lf\n",error_viterbi[0],(double)error_viterbi[0]/count);
		fprintf(fber,"viterbi M_ary error x : %d\tBER : %.12lf\n",error_viterbi_M_ary[0],(double)error_viterbi_M_ary[0]/count_M_ary);
		fprintf(fber_outSOVA, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_viterbi[1] / count);
		fprintf(fber_outSOVA_Mary, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_viterbi_M_ary[1] / count_M_ary);
		if (Encoder != 0)
		{
			fprintf(fber,"modulation error x : %d\tBER : %.12lf\n",error_demodul[0],(double)error_demodul[0]/count);
			fprintf(fber,"modulation M_ary error x : %d\tBER : %.12lf\n",error_demodul_M_ary[0],(double)error_demodul_M_ary[0]/count_modul_M_ary);
			fprintf(fber_outDEM, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_demodul[1] / count_modul);
			fprintf(fber_outDEM_Mary, "%d\t%.12lf\n", int(SNR_DB / SNR_INTERVAL), (double)error_demodul_M_ary[1] / count_modul_M_ary);
		}
		if (DRAWING_BER_CURVE == 0)
		{
			fprintf(BRIEF_RESULT, "%s_channel(%d)=%.12lf;\n", s_encoder[Encoder], int(SNR_DB / SNR_INTERVAL), (double)error / count);
			fprintf(BRIEF_RESULT, "%s_channel_M_ary(%d)=%.12lf;\n", s_encoder[Encoder], int(SNR_DB / SNR_INTERVAL), (double)error_M_ary / count_M_ary);
			fprintf(BRIEF_RESULT, "%%data(x)\n");
			//fprintf(BRIEF_RESULT,"%s_%s(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi[0]/count);
			//fprintf(BRIEF_RESULT,"%s_%s_M_ary(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi_M_ary[0]/count_M_ary);
			fprintf(BRIEF_RESULT, "%s_%s(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_viterbi[1] / count);
			fprintf(BRIEF_RESULT, "%s_%s_M_ary(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_viterbi_M_ary[1] / count_M_ary);
			if (Encoder != 0)
			{
				//fprintf(BRIEF_RESULT,"%s_%s_demodul(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul[0]/count_modul);
				//fprintf(BRIEF_RESULT,"%s_%s_demodul_M_ary(%d,:)=%.12lf;\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul_M_ary[0]/count_modul_M_ary);
				fprintf(BRIEF_RESULT, "%s_%s_demodul(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_demodul[1] / count_modul);
				fprintf(BRIEF_RESULT, "%s_%s_demodul_M_ary(%d,:)=%.12lf;\n", s_encoder[Encoder], s_detecter[DETECTOR], int(SNR_DB / SNR_INTERVAL), (double)error_demodul_M_ary[1] / count_modul_M_ary);
			}
		}
		
		if(DRAWING_BER_CURVE==1)
		{
			fprintf(BRIEF_RESULT,"%s_channel(%d)=%.12lf;\n",s_encoder[Encoder],int(SNR_DB/SNR_INTERVAL),(double)error/count);
			fprintf(BRIEF_RESULT,"%s_channel_M_ary(%d)=%.12lf;\n",s_encoder[Encoder],int(SNR_DB/SNR_INTERVAL),(double)error_M_ary/count_M_ary);
			fprintf(BRIEF_RESULT,"%%data(x,y)\n");
			fprintf(BRIEF_RESULT,"%s_%s(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi[0]/count,(double)error_viterbi[1]/count);
			fprintf(BRIEF_RESULT,"%s_%s_M_ary(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_viterbi_M_ary[0]/count_M_ary,(double)error_viterbi_M_ary[1]/count_M_ary);
			if(Encoder != 0)
			{
				fprintf(BRIEF_RESULT,"%s_%s_demodul(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul[0]/count_modul,(double)error_demodul[1]/count_modul);
				fprintf(BRIEF_RESULT,"%s_%s_demodul_M_ary(%d,:)=[%.12lf %.12lf];\n",s_encoder[Encoder],s_detecter[DETECTOR],int(SNR_DB/SNR_INTERVAL),(double)error_demodul_M_ary[0]/count_modul_M_ary,(double)error_demodul_M_ary[1]/count_modul_M_ary);
			}
		}
		
	}


	c_time = clock() - c_time; // c_time Check End 
	//fprintf(fber, "%d½Ă %dºĐ %dĂÊÀÔ´Ï´Ù. \n", int(c_time / CLOCKS_PER_SEC) / 3600, (int(c_time / CLOCKS_PER_SEC) % 3600) / 60, int(c_time / CLOCKS_PER_SEC) % 60);
	printf("%dhour %dmin %dsec. \n", int(c_time / CLOCKS_PER_SEC) / 3600, (int(c_time / CLOCKS_PER_SEC) % 3600) / 60, int(c_time / CLOCKS_PER_SEC) % 60);

	getchar();

	fcloseall();

	for (i = 0; i < SELL; i++)
		free(h[i]);
	free(h);

	for (i = 0; i < Target_Order; i++)
		free(Target_Coef[i]);
	free(Target_Coef);

	for (i = 0; i < Page_Size; i++)
		free(inLDPC[i]), free(EI[i]), free(outLDPC[i]);
	free(inLDPC), free(EI), free(outLDPC);

	for (i = 0; i < pow((double)M_ary, Target_Order - 1); i++)
		free(trellis[i]);
	free(trellis);

	for (i = 0; i < Equalizer_Length; i++)
		free(Equalizer_Coef_x[i]), free(Equalizer_Coef_y[i]);
	free(Equalizer_Coef_x), free(Equalizer_Coef_y);

	for (i = 0; i < Page_Size + (Equalizer_Length / 2) * 2; i++)
		free(inCH[i]), free(outCH[i]);
	free(inCH), free(outCH);

	for (i = 0; i < Page_Size; i++)
		free(outEQ_x[i]), free(outEQ_y[i]), free(n_outViterbi[i]), free(d_outViterbi[i]), free(d_outViterbi_x[i]), free(d_outViterbi_y[i]);
	free(d_outViterbi), free(n_outViterbi), free(outEQ_x), free(outEQ_y);
	free(d_outViterbi_x), free(d_outViterbi_y);


	free(input_park_2_3_4ary), free(output_park_2_3_4ary);

	for (i = 0; i < M_ary - 2; i++)
		free(trellis_A[i]);
	free(trellis_A);

	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			free(trellis_modul[i][j]);
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		free(trellis_modul[i]);
	free(trellis_modul);

	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		for (j = 0; j < (int)pow((double)M_ary, 2); j++)
			free(trellis_modul_prof[i][j]);
	for (i = 0; i < (int)pow((double)M_ary, 2); i++)
		free(trellis_modul_prof[i]);
	free(trellis_modul_prof);

	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			for (k = 0; k < 2; k++)
				free(trellis_modul_prof_method0[i][j][k]);
	for (i = 0; i < total_states_method0; i++)
		for (j = 0; j < total_states_method0; j++)
			free(trellis_modul_prof_method0[i][j]);
	for (i = 0; i < total_states_method0; i++)
		free(trellis_modul_prof_method0[i]);
	free(trellis_modul_prof_method0);




}


void HoloChannel(double** h, double SigmaB, double Mis_alignment_sigma, int SELL, int Unit_Sell)
{
	int i, j, l, m;

	double mis_length = gaussian(0.0, Mis_alignment_sigma, &k_noise);
	double mis_radian = gaussian(0.0, Mis_alignment_sigma, &k_noise);
	double mis_x = mis_length * cos(mis_radian);
	double mis_y = mis_length * sin(mis_radian);

	if (Mis_const == 1)
		mis_x = Mis_alignment_x, mis_y = Mis_alignment_y;;


	double** temp_h = (double**)malloc(SELL * Unit_Sell * sizeof(double*));
	for (i = 0; i < SELL * Unit_Sell; i++)
		temp_h[i] = (double*)calloc(SELL * Unit_Sell, sizeof(double));

	for (i = -SELL / 2; i <= SELL / 2; i++)
	{
		for (j = -SELL / 2; j <= SELL / 2; j++)
		{
			for (l = -Unit_Sell / 2; l <= Unit_Sell / 2; l++)
			{
				for (m = -Unit_Sell / 2; m <= Unit_Sell / 2; m++)
				{
					if ((i + mis_x) * Unit_Sell + l == 0 && (j + mis_y) * Unit_Sell + m == 0)
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB);
					else if ((i + mis_x) * Unit_Sell + l == 0 && (j + mis_y) * Unit_Sell + m != 0)
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB) * pow((sin(PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB) / (PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB)), 2);
					else if ((i + mis_x) * Unit_Sell + l != 0 && (j + mis_y) * Unit_Sell + m == 0)
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB) * pow((sin(PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB) / (PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB)), 2);
					else
						temp_h[i * Unit_Sell + l + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)][j * Unit_Sell + m + (Unit_Sell * (SELL / 2) + Unit_Sell / 2)] = 1 / (SigmaB * SigmaB) * pow((sin(PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB) / (PI * (i + l / float(Unit_Sell) + mis_x) / SigmaB)) * (sin(PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB) / (PI * (j + m / float(Unit_Sell) + mis_y) / SigmaB)), 2);
				}
			}
		}
	}
	for (i = 0; i < SELL; i++)
		for (j = 0; j < SELL; j++)
		{
			h[i][j] = 0;
			for (l = 0; l < Unit_Sell; l++)
				for (m = 0; m < Unit_Sell; m++)
					h[i][j] += temp_h[i * Unit_Sell + l][j * Unit_Sell + m];
		}

	for (i = 0; i < SELL * Unit_Sell; i++)
		free(temp_h[i]);
	free(temp_h);
}
void PR_Target_set(double** Target_Coef, int target_order)											//// Ă¤³ÎÀÇ ¸đ¾ç¿¡ µû¶ó PR TargetÀÇ coefficient¸¦ Á¤ÇÑ´Ù.
{
	if (target_order == 3)
	{
		Target_Coef[0][0] = 0, Target_Coef[0][1] = 1, Target_Coef[0][2] = 0;
		Target_Coef[1][0] = 1, Target_Coef[1][1] = 3, Target_Coef[1][2] = 1;
		//		if(DETECTOR==1)
		//			Target_Coef[1][0]=1, Target_Coef[1][1]=3, Target_Coef[1][2]=1;
		//		else if(DETECTOR==5)
		//			Target_Coef[1][0]=1, Target_Coef[1][1]=4, Target_Coef[1][2]=1;
		//		else
		//			Target_Coef[1][0]=1, Target_Coef[1][1]=4, Target_Coef[1][2]=1;
		//		Target_Coef[1][0]=1, Target_Coef[1][1]=int(Mis_alignment_dB), Target_Coef[1][2]=1;
		Target_Coef[2][0] = 0, Target_Coef[2][1] = 1, Target_Coef[2][2] = 0;
	}
	else if (target_order == 5)
	{
		Target_Coef[0][0] = 0, Target_Coef[0][1] = 0, Target_Coef[0][2] = 1, Target_Coef[0][3] = 0, Target_Coef[0][4] = 0;
		Target_Coef[1][0] = 0, Target_Coef[1][1] = 1, Target_Coef[1][2] = 2, Target_Coef[1][3] = 1, Target_Coef[1][4] = 0;
		Target_Coef[2][0] = 1, Target_Coef[2][1] = 2, Target_Coef[2][2] = 8, Target_Coef[2][3] = 2, Target_Coef[2][4] = 1;
		Target_Coef[3][0] = 0, Target_Coef[3][1] = 1, Target_Coef[3][2] = 2, Target_Coef[3][3] = 1, Target_Coef[3][4] = 0;
		Target_Coef[4][0] = 0, Target_Coef[4][1] = 0, Target_Coef[4][2] = 1, Target_Coef[4][3] = 0, Target_Coef[4][4] = 0;
	}
	else if (target_order == 7)
	{
		Target_Coef[0][0] = 0, Target_Coef[0][1] = 0, Target_Coef[0][2] = 0, Target_Coef[0][3] = 0, Target_Coef[0][4] = 0, Target_Coef[0][5] = 0, Target_Coef[0][6] = 0;
		Target_Coef[1][0] = 0, Target_Coef[1][1] = 0, Target_Coef[1][2] = 0, Target_Coef[1][3] = 0, Target_Coef[1][4] = 0, Target_Coef[1][5] = 0, Target_Coef[1][6] = 0;
		Target_Coef[2][0] = 0, Target_Coef[2][1] = 0, Target_Coef[2][2] = 0, Target_Coef[2][3] = 0, Target_Coef[2][4] = 0, Target_Coef[2][5] = 0, Target_Coef[2][6] = 0;
		Target_Coef[3][0] = 0, Target_Coef[3][1] = 0, Target_Coef[3][2] = 0, Target_Coef[3][3] = 0, Target_Coef[3][4] = 0, Target_Coef[3][5] = 0, Target_Coef[3][6] = 0;
		Target_Coef[4][0] = 0, Target_Coef[4][1] = 0, Target_Coef[4][2] = 0, Target_Coef[4][3] = 0, Target_Coef[4][4] = 0, Target_Coef[4][5] = 0, Target_Coef[4][6] = 0;
		Target_Coef[5][0] = 0, Target_Coef[5][1] = 0, Target_Coef[5][2] = 0, Target_Coef[5][3] = 0, Target_Coef[5][4] = 0, Target_Coef[5][5] = 0, Target_Coef[5][6] = 0;
		Target_Coef[6][0] = 0, Target_Coef[6][1] = 0, Target_Coef[6][2] = 0, Target_Coef[6][3] = 0, Target_Coef[6][4] = 0, Target_Coef[6][5] = 0, Target_Coef[6][6] = 0;
	}
}
void PRML_Convolution(double** trellis, double* Target_Coef, int Target_Order)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// trellis[i][j]	i=state, j= 0,1 => reference_bit(0 up, 1 down) 2,3 => reference_value(2 up, 3 down) ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Target_Order--;
	int i, j, k;					//// for¹® µ¹¸±¶§ ¾²À̀´Â ÀÓ½Ă º¯¼öµé
	int total_states = (int)pow((double)M_ary, Target_Order);


	int* temp_trel_a = (int*)calloc(Target_Order, sizeof(int));

	//// setting reference_bit ////
	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < M_ary; j++)
			trellis[i][j] = i % M_ary;
	}
	///////////////////////////////

	for (i = 0; i < Target_Order; i++)
		temp_trel_a[i] = 0;										//ĂÊ±âÈ­ÇØÁÖ´Â °Í 0
	for (i = 0; i < total_states; i++)
	{
		for (k = 0; k < Target_Order; k++)
			temp_trel_a[k] = (i / (int)pow((double)M_ary, k)) % M_ary;
		for (k = 0; k < M_ary; k++)
		{
			trellis[i][M_ary + k] = 0;
			for (j = 0; j < Target_Order; j++)
				trellis[i][M_ary + k] += Target_Coef[j] * temp_trel_a[j];
			trellis[i][M_ary + k] += (Target_Coef[Target_Order] * k);
			trellis[i][M_ary + k] /= M_ary - 1;
		}
	}


	free(temp_trel_a);
}
void Hologram_ChannelOut(double** outCH, double** h, int** inCH, int Page_Size, int SELL, double sigma)
{
	int i, j, x, y;
	double out;

	double sigma_AWGN = sqrt(sigma);

	for (x = 0; x < Page_Size; x++)
	{
		for (y = 0; y < Page_Size; y++)
		{
			out = 0;
			for (i = -SELL / 2; i <= SELL / 2; i++)
			{
				for (j = -SELL / 2; j <= SELL / 2; j++)
				{
					if (x + i >= 0 && y + j >= 0 && x + i < Page_Size && y + j < Page_Size)
					{
						out += inCH[x + i][y + j] * h[i + SELL / 2][j + SELL / 2];
					}
				}
			}
			if (noise_free == 0)
				outCH[x][y] = out / Unit_Sell / (M_ary - 1);
			else if (noise_free == 1)
				outCH[x][y] = out / Unit_Sell / (M_ary - 1) + gaussian(0.0, sigma_AWGN, &k_noise);
		}
	}
}
void Hologram_EqualizerOut_x(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order)
{
	int i, j, l, m;

	double Target_Value;

	//FILE *ERROR=fopen("eqout_error.dat","a");
	for (i = 0; i < Page_Size; i++)
	{
		for (j = 0; j < Page_Size; j++)
		{
			outEQ[i][j] = 0;
			for (l = 0; l < Equalizer_Length; l++)
				//l=Equalizer_Length/2;
				for (m = 0; m < Equalizer_Length; m++)
					outEQ[i][j] += inEQ[i + l][j + m] * eq_coef[l][m];
			Target_Value = 0;
			if (Training == 0 && i > Equalizer_Length / 2 && j > Equalizer_Length / 2)
			{
				l = 0;
				//for(l=-Target_Order/2 ; l<=Target_Order/2 ; l++)
				for (m = -Target_Order / 2; m <= Target_Order / 2; m++)
					//						Target_Value += (inCH[i+Equalizer_Length/2+l][j+Equalizer_Length/2+m]==0 ? -1 : 1)*Target_Coef[l+Target_Order/2][m+Target_Order/2];
					Target_Value += (double)inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] / (M_ary - 1) * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];
				//Target_Value += (double)inCH[i+Equalizer_Length/2+l][j+Equalizer_Length/2+m]*Target_Coef[l+Target_Order/2][m+Target_Order/2];

				for (l = 0; l < Equalizer_Length; l++)
					//l=Equalizer_Length/2;
					for (m = 0; m < Equalizer_Length; m++)
						eq_coef[l][m] += inEQ[i + l][j + m] * (Target_Value - outEQ[i][j]) * DELTA;
			}
			//if(Training == 0)
			//	fprintf(ERROR,"%lf\n",Target_Value-outEQ[i][j]);			
		}
	}
	//fclose(ERROR);
}
void Hologram_EqualizerOut_y(double** outEQ, double** eq_coef, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int Training, int** inCH, double** Target_Coef, int Target_Order)
{
	int i, j, l, m;

	double Target_Value;

	//	FILE *ERROR=fopen("eqout_error.dat","w");
	for (j = 0; j < Page_Size; j++)
	{
		for (i = 0; i < Page_Size; i++)
		{
			outEQ[i][j] = 0;
			for (l = 0; l < Equalizer_Length; l++)
				//l=Equalizer_Length/2;
				for (m = 0; m < Equalizer_Length; m++)
					outEQ[i][j] += inEQ[i + l][j + m] * eq_coef[l][m];
			Target_Value = 0;
			if (Training == 0 && i > Equalizer_Length / 2 && j > Equalizer_Length / 2)
			{
				m = 0;
				for (l = -Target_Order / 2; l <= Target_Order / 2; l++)
					//	for(m=-Target_Order/2 ; m<=Target_Order/2 ; m++)
	//						Target_Value += (inCH[i+Equalizer_Length/2+l][j+Equalizer_Length/2+m]==0 ? -1 : 1)*Target_Coef[l+Target_Order/2][m+Target_Order/2];
					Target_Value += (double)inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] / (M_ary - 1) * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];

				for (l = 0; l < Equalizer_Length; l++)
					//l=Equalizer_Length/2;
					for (m = 0; m < Equalizer_Length; m++)
						eq_coef[l][m] += inEQ[i + l][j + m] * (Target_Value - outEQ[i][j]) * DELTA;
			}
			//			if(Training == 1)
			//				fprintf(ERROR,"%lf\n",Target_Value-outEQ[i][j]);			
		}
	}
}
void NPML_Pridictor_x(double* Pridictor_Coef_x, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length)
{
	int i, j, l, m;

	double Target_Value;
	double prierror, pout;

	//	int count=0;
	//	FILE *OUT, *OUT1;

	//	OUT=fopen("npin.dat","a");
	//	OUT1=fopen("nperror.dat","a");

	double* Pridictor_Memory = (double*)calloc((Pridictor_Length + 1), sizeof(double));

	for (i = 0; i < Page_Size; i++)
	{
		for (j = 0; j < Page_Size; j++)
		{
			pout = 0.0;
			for (l = 0; l < Pridictor_Length + 1; l++)
				pout += Pridictor_Memory[l] * Pridictor_Coef_x[l];

			Target_Value = 0;
			for (l = -Target_Order / 2; l < Target_Order / 2; l++)
				for (m = -Target_Order / 2; m < Target_Order / 2; m++)
					Target_Value += inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];

			prierror = inEQ[i][j] - Target_Value - pout;

			//			if(count>100000 && count<500000)
			//			{
			//				fprintf(OUT,"%lf\n",inEQ[i][j] - Target_Value);
			//				fprintf(OUT1,"%lf\n",prierror);
			//			}
			//			count++;


			for (l = 0; l < Pridictor_Length + 1; l++)
				Pridictor_Coef_x[l] += DELTA * prierror * Pridictor_Memory[l];

			///////////////////	 Pridictor_Memory ¾¾ÇÁÆ®   ////////////////////////
			for (l = Pridictor_Length; l > 0; l--)
				Pridictor_Memory[l] = Pridictor_Memory[l - 1];
			Pridictor_Memory[0] = inEQ[i][j] - Target_Value;
		}
	}
	free(Pridictor_Memory);
}
void NPML_Pridictor_y(double* Pridictor_Coef_y, double** inEQ, int Page_Size, int SELL, int Equalizer_Length, int** inCH, double** Target_Coef, int Target_Order, int Pridictor_Length)
{
	int i, j, l, m;

	double Target_Value;
	double prierror, pout;

	double* Pridictor_Memory = (double*)calloc((Pridictor_Length + 1), sizeof(double));

	for (j = 0; j < Page_Size; j++)
	{
		for (i = 0; i < Page_Size; i++)
		{
			pout = 0.0;
			for (l = 0; l < Pridictor_Length + 1; l++)
				pout += Pridictor_Memory[l] * Pridictor_Coef_y[l];

			Target_Value = 0;
			for (l = -Target_Order / 2; l < Target_Order / 2; l++)
				for (m = -Target_Order / 2; m < Target_Order / 2; m++)
					Target_Value += inCH[i + Equalizer_Length / 2 + l][j + Equalizer_Length / 2 + m] * Target_Coef[l + Target_Order / 2][m + Target_Order / 2];

			prierror = inEQ[i][j] - Target_Value - pout;

			for (l = 0; l < Pridictor_Length + 1; l++)
				Pridictor_Coef_y[l] += DELTA * prierror * Pridictor_Memory[l];

			///////////////////	 Pridictor_Memory ¾¾ÇÁÆ®   ////////////////////////
			for (l = Pridictor_Length; l > 0; l--)
				Pridictor_Memory[l] = Pridictor_Memory[l - 1];
			Pridictor_Memory[0] = inEQ[i][j] - Target_Value;
		}
	}
	free(Pridictor_Memory);
}
void Hologram_SOVA_x(double** viterbi_out, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth)
{
	Target_Order--;
	int i, j, k, x, y;				//// for¹® µ¹¸±¶§ ¾²À̀´Â ÀÓ½Ă º¯¼öµé
	int Viterbi_temp;
	int vout = 0;
	double pri_error[M_ary];
	int over;
	double PM[M_ary];
	double alpha = 1.0;
	int error_viterbi = 0;
	int scale = 30;
	int iteration = 3;
	int total_states = (int)pow((double)M_ary, Target_Order);

	int ggg = 0;


	double eqout_total = 0, eqout_average;
	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			eqout_total += eq_out[i][j];
	eqout_average = eqout_total / (Page_Size * Page_Size);

	double** differ = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		differ[i] = (double*)calloc(M_ary, sizeof(double));

	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			viterbi_out[i][j] = 10000;

	unsigned long int*** Viterbi_State = (unsigned long int***)malloc(total_states * sizeof(unsigned long int**));
	for (i = 0; i < total_states; i++)
		Viterbi_State[i] = (unsigned long int**)malloc(depth * sizeof(unsigned long int*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_State[i][j] = (unsigned long int*)calloc(2, sizeof(unsigned long int));

	double** Viterbi_Path = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		Viterbi_Path[i] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_whereState = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_whereState[i] = (double**)malloc(depth * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_whereState[i][j] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_error_Memory = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_error_Memory[i] = (double**)malloc((Pridictor_Length + 2) * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			Viterbi_error_Memory[i][j] = (double*)calloc(2, sizeof(double));

	for (x = 0; x < Page_Size; x++)
	{
		for (i = 1; i < total_states; i++)
			Viterbi_Path[i][1] = 10000, Viterbi_State[i][0][0] = 0, Viterbi_State[i][0][1] = 0, Viterbi_State[i][1][0] = 0, Viterbi_State[i][1][1] = 0;
		Viterbi_Path[0][1] = 0, Viterbi_State[0][0][0] = 0, Viterbi_State[0][0][1] = 0, Viterbi_State[0][1][0] = 0, Viterbi_State[0][1][1] = 0;

		for (i = 0; i < total_states; i++)
			for (j = 0; j < Pridictor_Length + 2; j++)
				Viterbi_error_Memory[i][j][0] = 0.0, Viterbi_error_Memory[i][j][1] = 0.0;

		for (y = 0; y < Page_Size + 1; y++)
		{

			ggg = 0;
			for (i = 0; i < total_states; i++)
			{
				// calculate the predicted error
				pri_error[0] = 0.0, pri_error[1] = 0.0;
				for (j = 1; j < Pridictor_Length + 2; j++)
				{
					for (k = 0; k < M_ary; k++)
						pri_error[k] = Viterbi_error_Memory[i / M_ary + k * M_ary][j][1] * Pridictor_Coef[j - 1];
				}
				// calculate each path metric of a states
				if (NoiseFilter == 0)	for (k = 0; k < M_ary; k++) pri_error[k] = 0.0;
				if (y < Page_Size + 1 && y>0)
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1] + pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}
				else
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1];// +pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}

				// find minimum path of each state				
				Viterbi_temp = 0;
				if (gg_pr_viterbi == 1)
				{
					for (k = 0; k < M_ary; k++)
					{
						if (k == (M_ary - 1))
						{
							if ((i > 3) && (i < 13))
							{
								if (PM[Viterbi_temp] > PM[k])
									Viterbi_temp = k;
							}
						}
						else
						{
							if (PM[Viterbi_temp] > PM[k])
								Viterbi_temp = k;
						}
					}
				}
				else
				{
					for (k = 0; k < M_ary; k++)
						if (PM[Viterbi_temp] > PM[k])
							Viterbi_temp = k;
				}


				// save differents between survivor path and competition path
				for (k = 0; k < M_ary; k++)
					differ[i][k] = PM[k] - PM[Viterbi_temp];

				// save the survivor path of each state
				Viterbi_Path[i][0] = PM[Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_State[i][j][0] = Viterbi_State[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_State[i][0][0] = (unsigned long)trellis[i][Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_whereState[i][j][0] = Viterbi_whereState[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_whereState[i][0][0] = Viterbi_temp;
				for (j = 1; j < Pridictor_Length + 2; j++)
					Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / M_ary + Viterbi_temp * M_ary][j][1];
				if (y < Page_Size + 1 && y>0)	Viterbi_error_Memory[i][0][0] = eq_out[x][y - 1] - trellis[i][M_ary + Viterbi_temp];
			}

			// trace back and decide viterbi out
			if (y >= (Target_Order))
			{
				Viterbi_temp = 0;
				if (gg_pr_viterbi == 1)
				{
					for (i = 0; i < total_states; i++)
					{
						if ((i != 3) || (i != 13))
						{
							if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
								Viterbi_temp = i;
						}
					}
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);

				}
				else
				{
					for (i = 0; i < total_states; i++)
						if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
							Viterbi_temp = i;
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
				}
			}

			// shift the memory of path and state and errors¸̃¸đ¸®¸¦ shift½ĂÅ°´Â ±¸Á¶
			for (i = 0; i < total_states; i++)
			{
				Viterbi_Path[i][1] = Viterbi_Path[i][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
				for (j = 0; j < Pridictor_Length + 1; j++)
					Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
			}

			// prevent the path metric overflow
			// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
			over = 0;
			for (i = 0; i < total_states; i++)
				if (Viterbi_Path[i][1] > 10000)	over++;

			if (over == total_states)
				for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
		}

		for (i = depth - 1; i >= 1; i--)
			viterbi_out[x][Page_Size - (i)] = int(Viterbi_State[Viterbi_temp][i][0]);

	}

	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			free(Viterbi_State[i][j]), free(Viterbi_whereState[i][j]);
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			free(Viterbi_error_Memory[i][j]);
	for (i = 0; i < total_states; i++)
		free(Viterbi_State[i]), free(Viterbi_Path[i]), free(Viterbi_whereState[i]), free(Viterbi_error_Memory[i]);
	free(Viterbi_State), free(Viterbi_Path), free(Viterbi_whereState), free(Viterbi_error_Memory);
	free(differ);

}
void Hologram_SOVA_y(double** viterbi_out, double** eq_out, double** EI, double** trellis, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth)
{
	Target_Order--;
	int i, j, k, x, y;				//// for¹® µ¹¸±¶§ ¾²À̀´Â ÀÓ½Ă º¯¼öµé
	int Viterbi_temp;
	int vout = 0;
	double pri_error[M_ary];
	int over;
	double PM[M_ary];
	double alpha = 1.0;
	int error_viterbi = 0;
	int scale = 30;
	int iteration = 3;
	int total_states = (int)pow((double)M_ary, Target_Order);

	double eqout_total = 0, eqout_average;
	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			eqout_total += eq_out[i][j];
	eqout_average = eqout_total / (Page_Size * Page_Size);

	double** differ = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		differ[i] = (double*)calloc(M_ary, sizeof(double));

	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			viterbi_out[i][j] = 10000;

	unsigned long int*** Viterbi_State = (unsigned long int***)malloc(total_states * sizeof(unsigned long int**));
	for (i = 0; i < total_states; i++)
		Viterbi_State[i] = (unsigned long int**)malloc(depth * sizeof(unsigned long int*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_State[i][j] = (unsigned long int*)calloc(2, sizeof(unsigned long int));

	double** Viterbi_Path = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		Viterbi_Path[i] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_whereState = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_whereState[i] = (double**)malloc(depth * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_whereState[i][j] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_error_Memory = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_error_Memory[i] = (double**)malloc((Pridictor_Length + 2) * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			Viterbi_error_Memory[i][j] = (double*)calloc(2, sizeof(double));

	for (y = 0; y < Page_Size; y++)
	{
		for (i = 1; i < total_states; i++)
			Viterbi_Path[i][1] = 10000, Viterbi_State[i][0][0] = 0, Viterbi_State[i][0][1] = 0, Viterbi_State[i][1][0] = 0, Viterbi_State[i][1][1] = 0;
		Viterbi_Path[0][1] = 0, Viterbi_State[0][0][0] = 0, Viterbi_State[0][0][1] = 0, Viterbi_State[0][1][0] = 0, Viterbi_State[0][1][1] = 0;

		for (i = 0; i < total_states; i++)
			for (j = 0; j < Pridictor_Length + 2; j++)
				Viterbi_error_Memory[i][j][0] = 0.0, Viterbi_error_Memory[i][j][1] = 0.0;

		for (x = 0; x < Page_Size + 1; x++)
		{
			for (i = 0; i < total_states; i++)
			{
				// calculate the predicted error
				pri_error[0] = 0.0, pri_error[1] = 0.0;
				for (j = 1; j < Pridictor_Length + 2; j++)
				{
					for (k = 0; k < M_ary; k++)
						pri_error[k] = Viterbi_error_Memory[i / M_ary + k * M_ary][j][1] * Pridictor_Coef[j - 1];
				}
				// calculate each path metric of a states
				if (NoiseFilter == 0)	for (k = 0; k < M_ary; k++) pri_error[k] = 0.0;
				if (x < Page_Size + 1 && x>0)
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1] + pow(eq_out[x - 1][y] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}
				else
				{
					for (k = 0; k < M_ary; k++)
						PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1];// +pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
				}
				// find minimum path of each state
				Viterbi_temp = 0;
				for (k = 0; k < M_ary; k++)
					if (PM[Viterbi_temp] > PM[k])	Viterbi_temp = k;

				// save differents between survivor path and competition path
				for (k = 0; k < M_ary; k++)
					differ[i][k] = PM[k] - PM[Viterbi_temp];

				// save the survivor path of each state
				Viterbi_Path[i][0] = PM[Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_State[i][j][0] = Viterbi_State[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_State[i][0][0] = (unsigned long)trellis[i][Viterbi_temp];
				for (j = 1; j < depth; j++)
					Viterbi_whereState[i][j][0] = Viterbi_whereState[i / M_ary + Viterbi_temp * M_ary][j][1];
				Viterbi_whereState[i][0][0] = Viterbi_temp;
				for (j = 1; j < Pridictor_Length + 2; j++)
					Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / M_ary + Viterbi_temp * M_ary][j][1];
				if (x < Page_Size + 1 && x>0)	Viterbi_error_Memory[i][0][0] = eq_out[x - 1][y] - trellis[i][M_ary + Viterbi_temp];
			}

			// trace back and decide viterbi out
			if (x >= (Target_Order))
			{
				Viterbi_temp = 0;
				for (i = 0; i < total_states; i++)
					if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])	Viterbi_temp = i;

				if (x >= (depth - 1))
					viterbi_out[x - (depth - 1)][y] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
			}

			// shift the memory of path and state and errors¸̃¸đ¸®¸¦ shift½ĂÅ°´Â ±¸Á¶
			for (i = 0; i < total_states; i++)
			{
				Viterbi_Path[i][1] = Viterbi_Path[i][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
				for (j = depth - 1; j > 0; j--)
					Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
				for (j = 0; j < Pridictor_Length + 1; j++)
					Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
			}

			// prevent the path metric overflow
			// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
			over = 0;
			for (i = 0; i < total_states; i++)
				if (Viterbi_Path[i][1] > 10000)	over++;

			if (over == total_states)
				for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
		}

		for (i = depth - 1; i >= 1; i--)
			viterbi_out[Page_Size - (i)][y] = int(Viterbi_State[Viterbi_temp][i][0]);
	}

	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			free(Viterbi_State[i][j]), free(Viterbi_whereState[i][j]);
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			free(Viterbi_error_Memory[i][j]);
	for (i = 0; i < total_states; i++)
		free(Viterbi_State[i]), free(Viterbi_Path[i]), free(Viterbi_whereState[i]), free(Viterbi_error_Memory[i]);
	free(Viterbi_State), free(Viterbi_Path), free(Viterbi_whereState), free(Viterbi_error_Memory);
	free(differ);
}

void PRML_Convolution_gg_grad_23(double** trellis, double* Target_Coef, int Target_Order, int in_size, int* input)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//// trellis[i][j]	i=state, j= 0,1 => reference_bit(0 up, 1 down) 2,3 => reference_value(2 up, 3 down) ////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Target_Order--;
	int i, j, k;					//// for¹® µ¹¸±¶§ ¾²À̀´Â ÀÓ½Ă º¯¼öµé
	int total_states = (int)pow((double)in_size, Target_Order);
	int ggg = 0;
	int* temp_trel_a = (int*)calloc(Target_Order, sizeof(int));


	//// setting reference_bit ////
	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < in_size; j++)
			trellis[i][j] = input[i % in_size];
	}
	///////////////////////////////

	for (i = 0; i < Target_Order; i++)
		temp_trel_a[i] = input[0];										//ĂÊ±âÈ­ÇØÁÖ´Â °Í 0
	for (i = 0; i < total_states; i++)
	{
		temp_trel_a[0] = input[i % in_size];
		temp_trel_a[1] = input[i / in_size];

		for (k = 0; k < in_size; k++)
		{
			trellis[i][in_size + k] = 0;
			for (j = 0; j < Target_Order; j++)
				trellis[i][in_size + k] += Target_Coef[j] * temp_trel_a[j];
			trellis[i][in_size + k] += (Target_Coef[Target_Order] * input[k]);
			trellis[i][in_size + k] /= M_ary - 1;

		}
	}

	free(temp_trel_a);

}
void PRML_Convolution_for_modul_gg_grad_23(int*** trellis_modul)
{
	int i, j, k;
	int total_states = (int)pow((double)M_ary, 2);
	int temp_piece;
	int ggg = 0;

	int odd_piece[16][3] = { {0,1,1},{0,1,2},{0,1,3},{0,2,1},{0,2,2},{0,2,3},{1,1,1},{1,1,2},
						  {1,1,3},{1,2,1},{1,2,2},{1,2,3},{2,1,1},{2,1,2},{2,1,3},{2,2,1} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			temp_piece = (i + j) % 16;
			for (k = 0; k < 3; k++)
				trellis_modul[i][j][k] = odd_piece[temp_piece][k];
			ggg = 0;
		}
	}
}
void PRML_Convolution_for_modul_gg_grad_23_prof(int*** trellis_modul)
{

	int i, j, k;
	int total_states = (int)pow((double)M_ary, 2);
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{0,0,2},{0,0,3},{0,1,1},{0,1,2},{0,1,3},{0,2,2},{0,2,3},
					  {0,3,3},{1,1,2},{1,1,3},{1,2,2},{1,2,3},{1,3,3},{2,2,3},{2,3,3} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			temp_piece = (i + j) % 16;
			for (k = 0; k < 3; k++)
				trellis_modul[i][j][k] = piece[temp_piece][k];
			ggg = 0;
		}
	}
}
void PRML_Convolution_for_modul_gg_dmin2(int*** trellis_modul)
{

	int i, j, k;
	int total_states = (int)pow((double)M_ary, 2);
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{1,3,3},{0,0,3},{1,3,1},{0,1,0},{1,2,2},{0,1,2},{1,2,0},{0,2,1},{1,1,3},{0,2,3},{1,1,1},{0,3,0},{1,0,2},{0,3,2},{1,0,0} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			temp_piece = (i + j) % 16;
			for (k = 0; k < 3; k++)
				trellis_modul[i][j][k] = piece[temp_piece][k];
			ggg = 0;
		}
	}
}
void PRML_Convolution_for_modul_gg_grad_23_prof_method0(int**** trellis_modul)
{

	int i, j, k, l;
	int total_states = ((int)pow((double)M_ary, 2)) / 2;
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{0,0,2},{0,0,3},{0,1,1},{0,1,2},{0,1,3},{0,2,2},{0,2,3},
					  {0,3,3},{1,1,2},{1,1,3},{1,2,2},{1,2,3},{1,3,3},{2,2,3},{2,3,3} };


	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			for (k = 0; k < 2; k++)
			{
				temp_piece = ((i * 2) + (j * 2) + k) % 16;
				ggg = 0;
				for (l = 0; l < 3; l++)
					trellis_modul[i][j][k][l] = piece[temp_piece][l];
			}
		}
	}
}
void PRML_Convolution_for_modul_gg_dmin2_method0(int**** trellis_modul)
{

	int i, j, k, l;
	int total_states = ((int)pow((double)M_ary, 2)) / 2;
	int temp_piece;
	int ggg = 0;

	int piece[16][3] = { {0,0,1},{1,3,3},{0,0,3},{1,3,1},{0,1,0},{1,2,2},{0,1,2},{1,2,0},{0,2,1},{1,1,3},{0,2,3},{1,1,1},{0,3,0},{1,0,2},{0,3,2},{1,0,0} };

	for (i = 0; i < total_states; i++)
	{
		for (j = 0; j < total_states; j++)
		{
			for (k = 0; k < 2; k++)
			{
				temp_piece = ((i * 2) + (j * 2) + k) % 16;
				ggg = 0;
				for (l = 0; l < 3; l++)
					trellis_modul[i][j][k][l] = piece[temp_piece][l];
			}
		}
	}
}

void Hologram_SOVA_x_for_gg_grad(double** viterbi_out, double** eq_out, double** EI, double** trellis, double** trellisA, double* Pridictor_Coef, int Pridictor_Length, int Target_Order, int Page_Size, int depth)
{
	Target_Order--;
	int i, j, k, x, y;				//// for¹® µ¹¸±¶§ ¾²À̀´Â ÀÓ½Ă º¯¼öµé
	int Viterbi_temp;
	int vout = 0;
	double pri_error[M_ary];
	int over;
	double PM[M_ary];
	double alpha = 1.0;
	int error_viterbi = 0;
	int scale = 30;
	int iteration = 3;
	int total_states = (int)pow((double)M_ary, Target_Order);
	int test_line = 0;
	int ggg = 0;


	double eqout_total = 0, eqout_average;
	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			eqout_total += eq_out[i][j];
	eqout_average = eqout_total / (Page_Size * Page_Size);

	double** differ = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		differ[i] = (double*)calloc(M_ary, sizeof(double));

	for (i = 0; i < Page_Size; i++)
		for (j = 0; j < Page_Size; j++)
			viterbi_out[i][j] = 10000;

	unsigned long int*** Viterbi_State = (unsigned long int***)malloc(total_states * sizeof(unsigned long int**));
	for (i = 0; i < total_states; i++)
		Viterbi_State[i] = (unsigned long int**)malloc(depth * sizeof(unsigned long int*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_State[i][j] = (unsigned long int*)calloc(2, sizeof(unsigned long int));

	double** Viterbi_Path = (double**)malloc(total_states * sizeof(double*));
	for (i = 0; i < total_states; i++)
		Viterbi_Path[i] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_whereState = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_whereState[i] = (double**)malloc(depth * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			Viterbi_whereState[i][j] = (double*)calloc(2, sizeof(double));

	double*** Viterbi_error_Memory = (double***)malloc(total_states * sizeof(double**));
	for (i = 0; i < total_states; i++)
		Viterbi_error_Memory[i] = (double**)malloc((Pridictor_Length + 2) * sizeof(double*));
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			Viterbi_error_Memory[i][j] = (double*)calloc(2, sizeof(double));

	for (x = 0; x < Page_Size; x++)
	{
		for (i = 1; i < total_states; i++)
			Viterbi_Path[i][1] = 10000, Viterbi_State[i][0][0] = 0, Viterbi_State[i][0][1] = 0, Viterbi_State[i][1][0] = 0, Viterbi_State[i][1][1] = 0;
		Viterbi_Path[0][1] = 0, Viterbi_State[0][0][0] = 0, Viterbi_State[0][0][1] = 0, Viterbi_State[0][1][0] = 0, Viterbi_State[0][1][1] = 0;

		for (i = 0; i < total_states; i++)
			for (j = 0; j < Pridictor_Length + 2; j++)
				Viterbi_error_Memory[i][j][0] = 0.0, Viterbi_error_Memory[i][j][1] = 0.0;


		test_line = x % 3;
		if (test_line == 0)
		{
			total_states = (int)pow((double)2, Target_Order);

			for (y = 0; y < Page_Size + 1; y++)
			{
				ggg = 0;
				for (i = 0; i < total_states; i++)
				{
					// calculate the predicted error
					pri_error[0] = 0.0, pri_error[1] = 0.0;
					for (j = 1; j < Pridictor_Length + 2; j++)
					{
						for (k = 0; k < 2; k++)
							pri_error[k] = Viterbi_error_Memory[i / 2 + k * 2][j][1] * Pridictor_Coef[j - 1];
					}
					// calculate each path metric of a states
					if (NoiseFilter == 0)	for (k = 0; k < 2; k++) pri_error[k] = 0.0;

					if (y < Page_Size + 1 && y>0)
					{
						for (k = 0; k < 2; k++)
							PM[k] = Viterbi_Path[i / 2 + k * 2][1] + pow(eq_out[x][y - 1] - trellisA[i][2 + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}
					else
					{
						for (k = 0; k < 2; k++)
							PM[k] = Viterbi_Path[i / 2 + k * 2][1];//+pow(eq_out[x][y-1]-trellis[i][M_ary+k]-pri_error[0],2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}

					// find minimum path of each state				
					Viterbi_temp = 0;
					for (k = 0; k < 2; k++)
						if (PM[Viterbi_temp] > PM[k])
							Viterbi_temp = k;

					// save differents between survivor path and competition path
					for (k = 0; k < 2; k++)
						differ[i][k] = PM[k] - PM[Viterbi_temp];

					// save the survivor path of each state
					Viterbi_Path[i][0] = PM[Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_State[i][j][0] = Viterbi_State[i / 2 + Viterbi_temp * 2][j][1];
					Viterbi_State[i][0][0] = (unsigned long)trellisA[i][Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_whereState[i][j][0] = Viterbi_whereState[i / 2 + Viterbi_temp * 2][j][1];
					Viterbi_whereState[i][0][0] = Viterbi_temp;
					for (j = 1; j < Pridictor_Length + 2; j++)
						Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / 2 + Viterbi_temp * 2][j][1];
					if (y < Page_Size + 1 && y>0)	Viterbi_error_Memory[i][0][0] = eq_out[x][y - 1] - trellisA[i][2 + Viterbi_temp];
				}

				// trace back and decide viterbi out
				if (y >= (Target_Order))
				{
					Viterbi_temp = 0;
					for (i = 0; i < total_states; i++)
						if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
							Viterbi_temp = i;
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
				}

				// shift the memory of path and state and errors¸̃¸đ¸®¸¦ shift½ĂÅ°´Â ±¸Á¶
				for (i = 0; i < total_states; i++)
				{
					Viterbi_Path[i][1] = Viterbi_Path[i][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
					for (j = 0; j < Pridictor_Length + 1; j++)
						Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
				}

				// prevent the path metric overflow
				// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
				over = 0;
				for (i = 0; i < total_states; i++)
					if (Viterbi_Path[i][1] > 10000)	over++;

				if (over == total_states)
					for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
			}
		}
		else
		{
			total_states = (int)pow((double)4, Target_Order);

			for (y = 0; y < Page_Size + 1; y++)
			{
				ggg = 0;
				for (i = 0; i < total_states; i++)
				{
					// calculate the predicted error
					pri_error[0] = 0.0, pri_error[1] = 0.0;
					for (j = 1; j < Pridictor_Length + 2; j++)
					{
						for (k = 0; k < M_ary; k++)
							pri_error[k] = Viterbi_error_Memory[i / M_ary + k * M_ary][j][1] * Pridictor_Coef[j - 1];
					}
					// calculate each path metric of a states
					if (NoiseFilter == 0)	for (k = 0; k < M_ary; k++) pri_error[k] = 0.0;

					if (y < Page_Size + 1 && y>0)
					{
						for (k = 0; k < M_ary; k++)
							PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1] + pow(eq_out[x][y - 1] - trellis[i][M_ary + k] - pri_error[0], 2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}
					else
					{
						for (k = 0; k < M_ary; k++)
							PM[k] = Viterbi_Path[i / M_ary + k * M_ary][1];//+pow(eq_out[x][y-1]-trellis[i][M_ary+k]-pri_error[0],2);// - (EI[x][y-1]/2)*(trellis[i][k]*2-1);
					}

					// find minimum path of each state		

					Viterbi_temp = 0;
					for (k = 0; k < M_ary; k++)
						if (PM[Viterbi_temp] > PM[k])
							Viterbi_temp = k;

					// save differents between survivor path and competition path
					for (k = 0; k < M_ary; k++)
						differ[i][k] = PM[k] - PM[Viterbi_temp];

					// save the survivor path of each state
					Viterbi_Path[i][0] = PM[Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_State[i][j][0] = Viterbi_State[i / M_ary + Viterbi_temp * M_ary][j][1];
					Viterbi_State[i][0][0] = (unsigned long)trellis[i][Viterbi_temp];
					for (j = 1; j < depth; j++)
						Viterbi_whereState[i][j][0] = Viterbi_whereState[i / M_ary + Viterbi_temp * M_ary][j][1];
					Viterbi_whereState[i][0][0] = Viterbi_temp;
					for (j = 1; j < Pridictor_Length + 2; j++)
						Viterbi_error_Memory[i][j][0] = Viterbi_error_Memory[i / M_ary + Viterbi_temp * M_ary][j][1];
					if (y < Page_Size + 1 && y>0)	Viterbi_error_Memory[i][0][0] = eq_out[x][y - 1] - trellis[i][M_ary + Viterbi_temp];
				}

				// trace back and decide viterbi out
				if (y >= (Target_Order))
				{
					Viterbi_temp = 0;
					for (i = 0; i < total_states; i++)
						if (Viterbi_Path[Viterbi_temp][0] > Viterbi_Path[i][0])
							Viterbi_temp = i;
					if (y >= (depth - 1))
						viterbi_out[x][y - (depth - 1)] = int(Viterbi_State[Viterbi_temp][depth - 1][0]);
				}


				// shift the memory of path and state and errors¸̃¸đ¸®¸¦ shift½ĂÅ°´Â ±¸Á¶
				for (i = 0; i < total_states; i++)
				{
					Viterbi_Path[i][1] = Viterbi_Path[i][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_State[i][j][1] = Viterbi_State[i][j - 1][0];
					for (j = depth - 1; j > 0; j--)
						Viterbi_whereState[i][j][1] = Viterbi_whereState[i][j - 1][0];
					for (j = 0; j < Pridictor_Length + 1; j++)
						Viterbi_error_Memory[i][j + 1][1] = Viterbi_error_Memory[i][j][0];
				}

				// prevent the path metric overflow
				// If all path metric is greater than 10000, all path metric is subtracted from the 10000.
				over = 0;
				for (i = 0; i < total_states; i++)
					if (Viterbi_Path[i][1] > 10000)	over++;


				if (over == total_states)
					for (i = 0; i < total_states; i++)	Viterbi_Path[i][1] -= 10000;
			}

		}

		for (i = depth - 1; i >= 1; i--)
			viterbi_out[x][Page_Size - (i)] = int(Viterbi_State[Viterbi_temp][i][0]);


	}
	total_states = (int)pow((double)4, Target_Order);

	for (i = 0; i < total_states; i++)
		for (j = 0; j < depth; j++)
			free(Viterbi_State[i][j]), free(Viterbi_whereState[i][j]);
	for (i = 0; i < total_states; i++)
		for (j = 0; j < Pridictor_Length + 2; j++)
			free(Viterbi_error_Memory[i][j]);
	for (i = 0; i < total_states; i++)
		free(Viterbi_State[i]), free(Viterbi_Path[i]), free(Viterbi_whereState[i]), free(Viterbi_error_Memory[i]);
	free(Viterbi_State), free(Viterbi_Path), free(Viterbi_whereState), free(Viterbi_error_Memory);
	free(differ);

}

