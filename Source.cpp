#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <new>
#include <ctime>
#include <random>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string>

using namespace std;
/* 
	data_points and *_val is for full dynamics, other is for equilibriation
*/
const int data_points = 30;
const float min_val = 1;
const float max_val = 3.5;
const int iterations = 5000;

/* 
	Filenames and paths for calculations output
*/
const string filename = "test/kawasaki_fracdim_bigbox_m0_32.txt"; /* Old file for equilibriation output */
const string filepath = "competing_data/64-eq_0m_t100000_t100_Kawa0_bbJKpc_p"; /* Current file (with folder) for full dynamics output*/
const string ending = ".txt";

struct compEquilib {
	int N;
	float T1, T2, p, m;/*
	float E[iterations] = { 0 };
	float M[iterations] = { 0 };
	float I[iterations] = { 0 };
	int P[iterations] = { -1 };*/
	bool m1, m2, uni;
};

/* 
	Below is a datastruct to store spatial diversity values and print to file the full dynamics
*/
struct dataIsing { 
	float T1[data_points] = { 0 };
	float T2[data_points] = { -1 }; /* This means that kawasaki works as a pumping device */
	float E[data_points] = { 0 };
	float M[data_points] = { 0 };
	float M_up[data_points] = { 0 };
	float M_down[data_points] = { 0 };
	float I[data_points] = { 0 };
	float I_95up[data_points] = { 0 };
	float I_95down[data_points] = { 0 };
	float F[data_points] = { 0 };
	float F_95up[data_points] = { 0 };
	float F_95down[data_points] = { 0 };
};

/* 
	Below is a protoype class for single dynamics implementation
	Classes Metropolis and Kawasaki inherit this ones variables and methods
	Class named Competing is based on this one and functions are sometimes duplicate - more descriptions in that file
*/
class Ising {
public:
	int N, t;
	float beta, initM;

	/* lattice array is created dinamically
	extensive values that are yet to be averaged like E and M can be calculated from av variable or written in these arrays
	index lists are required for median and percentile calculations
	*/

	int* Energy;
	int* Mag;

	float* Index;
	float* I_95up;
	float* I_95down;

	float* FracDim;
	float* F_95up;
	float* F_95down;

	int** array;

	vector<int> up_x;
	vector<int> up_y;
	vector<int> down_x;
	vector<int> down_y;


	Ising(int n, float kT, float m, int iter, bool modelType, bool uniform) {

		N = n; /* Size of array */
		beta = 1 / kT; /* Temperature parameter Beta */
		initM = m; /* initial magnetization, remains const for Kawasaki */

		t = iter; /* Number of iterations, also check that number is even */
		if (t % 2 != 0) {
			t++;
		}

		/* memory allocation for lattice array */
		array = new int* [N];
		for (int i = 0; i < N; i++) {
			array[i] = new int[N];
		}

		/* memory allocation for energy and mag and index*/
		Energy = new int[t];
		Mag = new int[t];

		Index = new float[t];
		I_95up = new float[t];
		I_95down = new float[t];

		FracDim = new float[t];
		F_95up = new float[t];
		F_95down = new float[t];

		/* fill lattice */
		if (modelType == true) { /* true by default is for Metropolis, false for Kawasaki */
			if (uniform) {
				initUniform(initConfig); /*should be only - 1 or 1, chosen randomly */
			}
			else {
				initRandom();
			}
		}
		else {
			initKawasaki(initM);
		}
	}

	~Ising() {
		for (int i = 0; i < N; i++) {
			delete[] array[i];
		}
		delete[] array;

		delete[] Energy;
		delete[] Mag;

		delete[] Index;
		delete[] I_95up;
		delete[] I_95down;

		delete[] FracDim;
		delete[] F_95down;
		delete[] F_95up;
	}

	/* Initializing functions for random, uniform and kawasaki configurations*/
	void initUniform(int initConfig) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				array[i][j] = initConfig; // should be -1 or 1, default is 1
			}
		}
	}

	void initRandom() {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				array[i][j] = 2 * (rand() % 2) - 1; // 50/50
			}
		}
	}

	void initKawasaki(float m) {
		/* Create the correct amount of up/down spins in lattice to correspond to input m
		m is a number from -1 to 1, meaning total intensive magnetization*/
		int* darr = new int[N * N];
		random_device rd;
		mt19937 g(rd());

		for (int i = 0; i < N * N; i++) {
			if (i < N * N * 0.5 * (m + 1)) {
				darr[i] = 1;
			}
			else {
				darr[i] = -1;
			}
		}
		/* Shuffle the counted up/down spins and write into 2d array*/
		shuffle(darr, darr + N * N, g);

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				array[i][j] = darr[i * N + j];
			}
		}

		delete[] darr;
	}

	/* Functions to calculate energy and order parameter (magnetization)
	Also utility function to calculate neighbour spins is here*/
	int calcEnergy() {
		int temp_Energy = 0;

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				int s = array[i][j];
				int pos[2] = { i, j };
				int n = stateNeighbours(pos);

				temp_Energy += -n * s;
			}
		}

		return temp_Energy / 2;
	}

	int calcMag() {
		int temp_Mag = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				temp_Mag += array[i][j];
			}
		}

		return temp_Mag;
	}

	int stateNeighbours(int pos[]) {
		int xdown = pos[0] - 1;
		if (xdown < 0) {
			xdown = N - 1;
		}
		int ydown = pos[1] - 1;
		if (ydown < 0) {
			ydown = N - 1;
		}
		return array[(pos[0] + 1) % N][pos[1]] +
			array[xdown % N][pos[1]] +
			array[pos[0]][(pos[1] + 1) % N] +
			array[pos[0]][ydown % N];
	}


	/* Following two functions are for working with up/down spin lists
	vectors are employed to dynamically store and move spin position information*/
	void initializeLists() {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (array[i][j] == 1) {
					up_x.push_back(i);
					up_y.push_back(j);
				}
				else {
					down_x.push_back(i);
					down_y.push_back(j);
				}
			}
		}
	}

	void exchangeSpins(int up, int down) { /* x1,y1 for up spin list, x2,y2 for down spin list*/
		swap(up_x[up], down_x[down]);
		swap(up_y[up], down_y[down]);
	}

	/* Functions for diversity index I calculation. Performs scaling and calculating */

	float trapezoid(int ptr_scales[], float ptr_vals[], int arrLen) { /* Trapezoid 1/2(xi-xj)(fxi + fxj), j=i+1*/
		float integral = 0;

		for (int i = 0; i < arrLen - 1; i++) {
			integral = integral + (ptr_scales[i + 1] - ptr_scales[i]) * (ptr_vals[i + 1] + ptr_vals[i]) / 2;
		}

		return integral;
	}

	float autoScale(int scale) { /* returns a value of std for resized lattice to index calculation */

		/* Create new resized lattice */
		int N_new = ceil(float(N) / float(scale));
		float** resizedArr = new float* [N_new];

		for (int i = 0; i < N_new; i++) {
			resizedArr[i] = new float[N_new];
		}

		/* Fill the new resized lattice with mean value from contracted compartments
		Also calculate the mean for the resized lattice to have it for std dev calc*/
		float sumMeanArr = 0;

		for (int x = 0; x < N_new; x++) {
			for (int y = 0; y < N_new; y++) {

				float resizedSum = 0;
				for (int i = 0; i < scale; i++) {
					for (int j = 0; j < scale; j++) {
						resizedSum += array[scale * x + i][scale * y + j];
					}
				}
				resizedArr[x][y] = resizedSum / powf(scale, 2);
				sumMeanArr += resizedArr[x][y];
			}
		}
		sumMeanArr = sumMeanArr / (powf(N_new, 2));

		/* Calculate the standard deviation */
		float std = 0;

		for (int x = 0; x < N_new; x++) {
			for (int y = 0; y < N_new; y++) {
				std += powf((resizedArr[x][y] - sumMeanArr), 2);
			}
		}

		std = sqrt(std / powf(N_new, 2)); // ne scale o N_new????

		for (int i = 0; i < N_new; i++) {
			delete[] resizedArr[i];
		}

		return std;
	}

	float calcIndex() {

		/* Get scales and diversity measurement values */
		int scales_arrLen = floor(log2(N));
		int* scales = new int[scales_arrLen];
		float* vals = new float[scales_arrLen];
		float* minimal_vals = new float[scales_arrLen];
		float* baseline_vals = new float[scales_arrLen];

		for (int i = 0; i < scales_arrLen; i++) {
			scales[i] = pow(2, i);
			vals[i] = autoScale(scales[i]);
			minimal_vals[i] = 0;
			baseline_vals[i] = float(1 / float(scales[i])); /* changed to be unit size, pow(x, 2), corresponds to area! */
		}
		minimal_vals[0] = 1;

		/* Make diversity measurement values relative to the first one, accord to the finest scale*/
		if (vals[0] != 0) {
			for (int i = scales_arrLen - 1; i > -1; i--) {
				vals[i] = vals[i] / vals[0];
			}
		}
		else {
			for (int i = 0; i < scales_arrLen; i++) {
				vals[i] = 0;
			}
			vals[0] = 1;
		}

		/* x axis, scales, have area dimension so need to be squared. previous scales used for easier purpose in autoscaling*/
		for (int i = 0; i < scales_arrLen; i++) {
			scales[i] = scales[i] * scales[i];
		}

		/* Calculate data, minimal area and baseline (random) trapezoid integration value*/
		float minimal_area = trapezoid(scales, minimal_vals, scales_arrLen);
		float baseline_area = trapezoid(scales, baseline_vals, scales_arrLen);
		float data_area = trapezoid(scales, vals, scales_arrLen);

		float index;
		if (data_area < baseline_area) {
			index = (data_area - minimal_area) / (baseline_area - minimal_area) - 1;
		}
		else {
			index = (data_area - baseline_area) / (scales[scales_arrLen - 1] - 1 - baseline_area);
		}

		delete[] scales;
		delete[] vals;
		delete[] minimal_vals;
		delete[] baseline_vals;

		return index;
	}

	/* The following methods are for fractal dimension of lattice configuration calculations */

	void makeLogspace(float* data_values, int points, float min_exp, float max_exp, float base = 2) {
		float diff = (max_exp - min_exp) / (points - 1);
		for (int i = 0; i < points; i++) {
			data_values[i] = powf(base, min_exp + diff * (i));
		}
	}

	void makeLinspace(float* data_values, int points, float min_val, float max_val) {
		float diff = (max_val - min_val) / (points - 1);
		for (int i = 0; i < points; i++) {
			data_values[i] = min_val + diff * i;
		}
	}

	int check_bin(float* row, float* col) {
		for (int k = ceil(row[0]); k <= floor(row[1]); k++) {
			for (int l = ceil(col[0]); l <= floor(col[1]); l++) {
				if (((k <= N-1) && (k >= 0)) && ((l <= N-1) && (l >= 0))) {
					if (array[k][l] > 0) {
						return 1;
					}
				}
			}
		}
		return 0;
	}

	int bin_count(float scale, float offset) {
		int count = 0;
		for (float i = -offset; i < N; i += scale) {
			for (float j = -offset; j < N; j += scale) {
				float rows[2] = { i, i + scale };
				float cols[2] = { j, j + scale };

				count += check_bin(rows, cols);
				//std::cout << i << " " << i + scale << " " << j << " " << j + scale << "   " << count << endl;
			}
		}
		return count;
	}

	void linearRegression(float* coeffs, float* x_vals_inv, float* y_vals, int n_vals) {
		/* simple linear regression https://www.bragitoff.com/2015/09/c-program-to-linear-fit-the-data-using-least-squares-method/ */
		float xsum = 0, x2sum = 0, ysum = 0, xysum = 0;

		for (int i = 0; i < n_vals; i++) {
			xsum = xsum + (x_vals_inv[i]);
			ysum = ysum + y_vals[i];
			x2sum = x2sum + pow((x_vals_inv[i]), 2);
			xysum = xysum + (x_vals_inv[i]) * y_vals[i];
		}

		coeffs[0] = (n_vals * xysum - xsum * ysum) / (n_vals * x2sum - xsum * xsum); /* slope */
		coeffs[1] = (x2sum * ysum - xsum * xysum) / (x2sum * n_vals - xsum * xsum); /* intercept */

	}

	float calcFracDim(int n_scales=20, int n_offsets=5, float min_exp=1, float max_exp=2) {
		max_exp = floor(log2(N / 2));
		float dim[2] = { 0 };
		float* scales = new float[n_scales]; /* input param */
		float* bin_vals = new float[n_scales];

		makeLogspace(scales, n_scales, min_exp, max_exp);

		for (int i = 0; i < n_scales; i++) {
			int* temp = new int[n_offsets];
			float* offsets = new float[n_offsets]; /* input param */

			makeLinspace(offsets, n_offsets, 0.01, scales[i] - 0.01);

			for (int j = 0; j < n_offsets; j++) {
				temp[j] = bin_count(scales[i], offsets[j]);
			}

			bin_vals[i] = *min_element(temp, temp + n_offsets);

			delete[] offsets;
			delete[] temp;

			bin_vals[i] = log(bin_vals[i]);
			scales[i] = log((1 / scales[i]));
		}

		linearRegression(dim, scales, bin_vals, n_scales);

		delete[] bin_vals;
		delete[] scales;

		return dim[0];
	}


private:

	int initConfig = 1;

	void setInitConfig(int a) {
		initConfig = a;
	}
};

class Metropolis : public Ising {
public:
	using Ising::Ising;

	float exponents[2] = { exp(-4 * beta), exp(-8 * beta) };
	/*
	float av_energy = 0;
	float av_mag = 0;*/

	void eqAvIter(float* ptr_energy, float* ptr_mag, float* ptr_index, float* ptr_indexUp, float* ptr_indexDown, float* ptr_frac, float* ptr_fracUp, float* ptr_fracDown) {

		float av_energy = 0;
		float av_mag = 0;
		float av_index = 0;

		/* equilibriating system for 2000 iterations */
		for (int i = 0; i < 2000; i++) {
			updateStep();
		}

		/* accumulating desired values for t iterations*/
		for (int i = 0; i < t; i++) {
			updateStep();
			av_energy += calcEnergy();
			av_mag += calcMag();
			//av_index += calcIndex(); /* FOR MEDIAN VALUE USE THE ALlOCATED MEMORY OF Index AND SORT IT. TAKE MIDDLE VALUE*/
			Index[i] = calcIndex();
			FracDim[i] = calcFracDim();
		}

		*ptr_energy = av_energy / (powf(N, 2) * t);
		*ptr_mag = av_mag / (powf(N, 2) * t);

		sort(Index, Index + t);
		*ptr_index = Index[t / 2];
		*ptr_indexUp = Index[t / 20];
		*ptr_indexDown = Index[19 * t / 20];

		sort(FracDim, FracDim + t);
		*ptr_frac = FracDim[t / 2];
		*ptr_fracUp = FracDim[t / 20];
		*ptr_fracDown = FracDim[19 * t / 20];
	}

	void updateStep() {

		/* random gen should be moved to the function with iter*/

		random_device rd; // obtain a random number from hardware
		mt19937 gen(rd()); // seed the generator			
		uniform_int_distribution<> distr(0, N - 1); // define the range for pos
		uniform_real_distribution<> distr1(0, 1); // define the range for rand 0-1

		int n_steps = powf(N, 2);
		for (int i = 0; i < n_steps; i++) {
			int select_pos[2] = { distr(gen), distr(gen) };

			int dE = 2 * array[select_pos[0]][select_pos[1]] * stateNeighbours(select_pos);

			if (dE <= 0) {
				array[select_pos[0]][select_pos[1]] *= -1;
			}
			else if (distr1(gen) < exponents[dE / 4 - 1]) {
				array[select_pos[0]][select_pos[1]] *= -1;
			}
		}
	}
};

class Kawasaki : public Ising {
public:
	using Ising::Ising;

	float exponents[4] = { exp(-4 * beta), exp(-8 * beta), exp(-12 * beta), exp(-16 * beta) };

	/* Following two functions are for dynamics implementation */
	void eqAvIter(float* ptr_energy, float* ptr_mag, float* ptr_index, float* ptr_indexUp, float* ptr_indexDown, float* ptr_frac, float* ptr_fracUp, float* ptr_fracDown) {


		float av_energy = 0;
		float av_mag = 0;
		float av_index = 0;

		/* Initialize lists */
		initializeLists();

		for (int i = 0; i < 2000; i++) {
			updateStep();
		}

		for (int i = 0; i < t; i++) {
			updateStep();
			av_energy += calcEnergy();
			av_mag += calcMag();
			Index[i] = calcIndex();
			FracDim[i] = calcFracDim();
		}

		*ptr_energy = av_energy / (powf(N, 2) * t);
		*ptr_mag = av_mag / (powf(N, 2) * t);

		sort(Index, Index + t);
		*ptr_index = Index[t / 2];
		*ptr_indexUp = Index[t / 20];
		*ptr_indexDown = Index[19 * t / 20];

		sort(FracDim, FracDim + t);
		*ptr_frac = FracDim[t / 2];
		*ptr_fracUp = FracDim[t / 20];
		*ptr_fracDown = FracDim[19 * t / 20];
	}

	void updateStep() {

		if ((up_x.size() == 0) || (down_x.size() == 0)) {
			cout << " One of the vectors is empty and therefore distr_up or distr_down is crashes\n";
			return;
		}
		random_device rd; // obtain a random number from hardware
		mt19937 gen(rd()); // seed the generator		
		uniform_int_distribution<> distr_up(0, up_x.size() - 1); // two ranges for up and down vectors
		uniform_int_distribution<> distr_down(0, down_x.size() - 1); // two ranges for up and down vectors
		uniform_real_distribution<> distr1(0, 1); // define the range for rand [0-1)

		int n_steps = powf(N, 2);
		for (int i = 0; i < n_steps; i++) {
			int pos[2] = { distr_up(gen), distr_down(gen) };
			int select_pos_up[2] = { up_x[pos[0]], up_y[pos[0]] };
			int select_pos_down[2] = { down_x[pos[1]], down_y[pos[1]] };

			int dEu = array[select_pos_up[0]][select_pos_up[1]] * stateNeighbours(select_pos_up) + array[select_pos_down[0]][select_pos_down[1]] * stateNeighbours(select_pos_down);

			array[select_pos_up[0]][select_pos_up[1]] *= -1;
			array[select_pos_down[0]][select_pos_down[1]] *= -1;

			int dEv = array[select_pos_up[0]][select_pos_up[1]] * stateNeighbours(select_pos_up) + array[select_pos_down[0]][select_pos_down[1]] * stateNeighbours(select_pos_down);

			if (((dEu - dEv) <= 0) || (distr1(gen) < exponents[(dEu - dEv) / 4 - 1])) {
				exchangeSpins(pos[0], pos[1]);
			}
			else {
				array[select_pos_up[0]][select_pos_up[1]] *= -1;
				array[select_pos_down[0]][select_pos_down[1]] *= -1;
			}

		}
	}

};

/*
	Main class for competing dynamics Ising object generation
	Class includes:
		dynamics implementation (single and competing)
		state parameter calculations
		diversity index calculation methods
		fractal dimension via box-counting method calculations
*/
class Competing{
public:
	
	/* lattice array is created dinamically
	extensive values that are yet to be averaged like E and M can be calculated from av variable or written in these arrays
	index lists are required for median and percentile calculations
	*/

	int* Energy;
	float* Mag;
	float* Index;
	float* FracDim;

	/*
		Lattice variable
	*/
	int** array;

	/*
		dynamic arrays to store locations of up/down spins
		something similar to "sparse" matrix entry storing
	*/
	vector<int> up_x;
	vector<int> up_y;
	vector<int> down_x;
	vector<int> down_y;

	int N, t1, t2;
	float beta1, beta2, dynamics_parameter, initM;
	bool lattice_param;	

	Competing(int n, float kT1, float kT2, float p, float m, int iter1, int iter2, bool uniform) {
		/*
			Constructor to initialize simulation parameters: 
				size of array (one dimension of 2D lattice)
				certain kT1 and kT2 temperatures 
				selected dynamic selection parameter p
				initial magnetization m (for kawasaki initial ordering)
				number of iterations
				lattice_param defines initial ordering (i.e. (anti)ferromagnetic or paramagnetic)
			Initialization is sufficient to start calculations
		*/
		N = n; /* Size of array */

		beta1 = 1 / kT1; /* Temperature parameter Beta */
		beta2 = 1 / kT2;

		initM = m; /* initial magnetization, remains const for Kawasaki */
		dynamics_parameter = p; /* probability for first type of dynamic (Metropolis) to occur.*/

		lattice_param = uniform;

		t1 = iter1; /* Number of iterations for equilibriation */
		t2 = iter2; /* Number of iterations for averaging */
		/*
		if (t % 2 != 0) {
			t++;
		}*/

		/* memory allocation for lattice array */
		array = new int* [N];
		for (int i = 0; i < N; i++) {
			array[i] = new int[N];
		}

		/* memory allocation for energy and mag and index*/
		Energy = new int[t2];
		Mag = new float[t2];
		Index = new float[t2];
		FracDim = new float[t2];

		/* fill lattice */
		if (lattice_param) {
			initUniform(1); /* should be only - 1 or 1 */
		}
		else {
			initKawasaki(initM); /* initializes a fixed magnetization lattice, by default random */
		}

		initializeLists();

	}

	~Competing() {
		for (int i = 0; i < N; i++) {
			delete[] array[i];
		}
		delete[] array;
		delete[] Energy;
		delete[] Mag;
		delete[] Index;
		delete[] FracDim;
	}

	/* 
		Initializing functions for random, uniform and Kawasaki configurations
		Kawasaki configuration has initial mean magnetization m which does not change
	*/

	void initUniform(int initConfig) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				array[i][j] = initConfig; // should be -1 or 1, default is 1
			}
		}
	}

	void initRandom() {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				array[i][j] = 2 * (rand() % 2) - 1; // 50/50
			}
		}
	}

	void initKawasaki(float m) {
		/* Create the correct amount of up/down spins in lattice to correspond to input m
		m is a number from -1 to 1, meaning total intensive magnetization*/
		int* darr = new int[N * N];
		random_device rd;
		mt19937 g(rd());

		for (int i = 0; i < N * N; i++) {
			if (i < N * N * 0.5 * (m + 1)) {
				darr[i] = 1;
			}
			else {
				darr[i] = -1;
			}
		}
		/* Shuffle the counted up/down spins and write into 2d array*/
		shuffle(darr, darr + N * N, g);

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				array[i][j] = darr[i * N + j];
			}
		}

		delete[] darr;
	}

	/* 
		Functions to calculate energy and order parameter (magnetization)
		Also utility function to calculate neighbour spins is here
	*/

	int calcEnergy() {
		int temp_Energy = 0;

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				int s = array[i][j];
				int pos[2] = { i, j };
				int n = stateNeighbours(pos);

				temp_Energy += -n * s;
			}
		}

		return temp_Energy / 2;
	}

	int calcMag() {
		int temp_Mag = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				temp_Mag += array[i][j];
			}
		}

		return temp_Mag;
	}

	int stateNeighbours(int pos[]) {
		int xdown = pos[0] - 1;
		if (xdown < 0) {
			xdown = N - 1;
		}
		int ydown = pos[1] - 1;
		if (ydown < 0) {
			ydown = N - 1;
		}
		return array[(pos[0] + 1) % N][pos[1]] +
			array[xdown % N][pos[1]] +
			array[pos[0]][(pos[1] + 1) % N] +
			array[pos[0]][ydown % N];
	}

	/* 
		Following three functions are for working with up/down spin lists
		vectors are employed to dynamically store and move spin position information
		vectors serve like entries of "sparse" matrix
	*/

	void initializeLists() {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (array[i][j] == 1) {
					up_x.push_back(i);
					up_y.push_back(j);
				}
				else {
					down_x.push_back(i);
					down_y.push_back(j);
				}
			}
		}
	}

	void exchangeSpins(int up, int down) { /* x1,y1 for up spin list, x2,y2 for down spin list*/
		swap(up_x[up], down_x[down]);
		swap(up_y[up], down_y[down]);
	}

	void transferSpin(int* pos) {

		if (array[pos[0]][pos[1]] == 1) { /* from down to up list */
			for (int i = 0; i < down_x.size(); i++) {
				if ((down_x[i] == pos[0]) && (down_y[i]) == pos[1]) {
					up_x.push_back(down_x[i]);
					up_y.push_back(down_y[i]);
					down_x.erase(down_x.begin() + i);
					down_y.erase(down_y.begin() + i);
					return;
				}
			}
		}
		else { /* from up list to down list */
			for (int i = 0; i < up_x.size(); i++) {
				if ((up_x[i] == pos[0]) && (up_y[i]) == pos[1]) {
					down_x.push_back(up_x[i]);
					down_y.push_back(up_y[i]);
					up_x.erase(up_x.begin() + i);
					up_y.erase(up_y.begin() + i);
					return;
				}
			}
		}
	}

	/* 
		Functions for diversity index I calculation
		Performs scaling and calculating
	*/

	float trapezoid(int ptr_scales[], float ptr_vals[], int arrLen) {
		/*
			Calculates area under the curve of relative scaling law.
			Uses Trapezoid integration rule: 1/2(xi-xj)(fxi + fxj), j=i+1
				ptr_scales are x's
				ptr_vals are y's
				arrLen is total number of values
		*/
		float integral = 0;

		for (int i = 0; i < arrLen - 1; i++) {
			integral = integral + (ptr_scales[i + 1] - ptr_scales[i]) * (ptr_vals[i + 1] + ptr_vals[i]) / 2;
		}

		return integral;
	}

	float autoScale(int scale) { /* returns a value of std for resized lattice to index calculation */
		/* 
			Creates new resized lattice after scaling the original lattice by size scale
			variable "scale" can also be thought of as a resolution multiplier (reducing resolution)
			transformation is similar to renormalization technique, except new values are mean of joined values
		*/
		int N_new = ceil(float(N) / float(scale));
		float** resizedArr = new float* [N_new];

		for (int i = 0; i < N_new; i++) {
			resizedArr[i] = new float[N_new];
		}

		/* Fill the new resized lattice with mean value from contracted compartments
		Also calculate the mean for the resized lattice to have it for std dev calc*/
		float sumMeanArr = 0;

		for (int x = 0; x < N_new; x++) {
			for (int y = 0; y < N_new; y++) {

				float resizedSum = 0;
				for (int i = 0; i < scale; i++) {
					for (int j = 0; j < scale; j++) {
						resizedSum += array[scale * x + i][scale * y + j];
					}
				}
				resizedArr[x][y] = resizedSum / powf(scale, 2);
				sumMeanArr += resizedArr[x][y];
			}
		}
		sumMeanArr = sumMeanArr / (powf(N_new, 2));

		/* 
			Calculating the standard deviation of spins in such lattice for scaling curve
			This standard deviation is later normalized to original lattice one (it is always smaller)
		*/
		float std = 0;

		for (int x = 0; x < N_new; x++) {
			for (int y = 0; y < N_new; y++) {
				std += powf((resizedArr[x][y] - sumMeanArr), 2);
			}
		}

		std = sqrt(std / powf(N_new, 2)); // ne scale o N_new????

		for (int i = 0; i < N_new; i++) {
			delete[] resizedArr[i];
		}

		return std;
	}

	float calcIndex() {
		
		/*
			Diversity index I calculation function. I values range in [-1, 1]
			Employs auto_scale function to get standard deviation value for scaled lattice
			For a number of scales_arrLen points a scaling curve is created
			trapezoid function is used to calculate area under it and form diversity index values in appropriate range
		*/

		int scales_arrLen = floor(log2(N));
		int* scales = new int[scales_arrLen];
		float* vals = new float[scales_arrLen];
		float* minimal_vals = new float[scales_arrLen];
		float* baseline_vals = new float[scales_arrLen];

		for (int i = 0; i < scales_arrLen; i++) {
			scales[i] = pow(2, i);
			vals[i] = autoScale(scales[i]);
			minimal_vals[i] = 0;
			baseline_vals[i] = float(1 / float(scales[i])); /* changed to be unit size, pow(x, 2), corresponds to area! */
		}
		minimal_vals[0] = 1;

		/* 
			Make diversity measurement values relative to the first one, accord to the finest resolution (original) scale
			This ensures the curve reproducability
		*/
		if (vals[0] != 0) {
			for (int i = scales_arrLen - 1; i > -1; i--) {
				vals[i] = vals[i] / vals[0];
			}
		}
		else {
			for (int i = 0; i < scales_arrLen; i++) {
				vals[i] = 0;
			}
			vals[0] = 1;
		}

		/* x axis, scales, have area dimension so need to be squared. previous scales used for easier purpose in autoscaling*/
		for (int i = 0; i < scales_arrLen; i++) {
			scales[i] = scales[i] * scales[i];
		}

		/* Calculate data, minimal area and baseline (random) trapezoid integration value*/
		float minimal_area = trapezoid(scales, minimal_vals, scales_arrLen);
		float baseline_area = trapezoid(scales, baseline_vals, scales_arrLen);
		float data_area = trapezoid(scales, vals, scales_arrLen);

		float index;
		if (data_area < baseline_area) {
			index = (data_area - minimal_area) / (baseline_area - minimal_area) - 1;
		}
		else {
			index = (data_area - baseline_area) / (scales[scales_arrLen - 1] - 1 - baseline_area);
		}

		delete[] scales;
		delete[] vals;
		delete[] minimal_vals;
		delete[] baseline_vals;

		return index;
	}

	/* 
		The following methods are for fractal dimension of lattice configuration calculations
		Box-counting method is employed which makes use of scaling boxes to be fitted on lattice
		Since the scaling used is always finite due to finitiness of lattice, it can be also termed finite-size scaling law
		Note: usually there are no single best fit (multifractality emerges) so this only takes best value for linear regression
	*/

	void makeLogspace(float* data_values, int points, float min_exp, float max_exp, float base = 2) {
		/*
			Makes logspace array, useful to have linear points in a log-log plot
		*/
		float diff = (max_exp - min_exp) / (points - 1);
		for (int i = 0; i < points; i++) {
			data_values[i] = powf(base, min_exp + diff * (i));
		}
	}

	void makeLinspace(float* data_values, int points, float min_val, float max_val) {
		/*
			Distributes values linspace in an array
		*/
		float diff = (max_val - min_val) / (points - 1);
		for (int i = 0; i < points; i++) {
			data_values[i] = min_val + diff * i;
		}
	}

	int check_bin(float* row, float* col, int a, float m) {
		/*
			Checks if the current box has a spin of desired value inside of it
			Since we dont need a full values of 2D histogram, just a bollean type True or False value - we can return 0 or 1
			Other return methods may be used (which are commented below):
				return 1 only if the spins in the box have similar mean magnetization
				if the mean magnetization is non zero/one
				etc.
			Such methods are better at grasping interfaces and solid bodies of lattice
		*/
		float sum = 0;
		float it = 0;
		for (int k = ceil(row[0]); k <= floor(row[1]); k++) {
			for (int l = ceil(col[0]); l <= floor(col[1]); l++) {
				if (((k <= N-1) && (k >= 0)) && ((l <= N-1) && (l >= 0))) {
					//it++;
					//sum += array[k][l];
					if (array[k][l] == a) {
						//sum += array[k][l];
						return 1;
					}
				}
			}
		}
		return 0;
		/*
		if (abs(sum)/it != 1) {
			return 1;
		}
		else {
			return 0;
		}
		*/
	}

	int bin_count(float scale, float offset, int a, float m) {
		/*
			Function to form the 2D histogram
			scale is the side length of square box
			offset is the value to be offset from initial array position (for count minimizing)
		*/
		int count = 0;
		for (float i = -offset; i < N; i += scale) {
			for (float j = -offset; j < N; j += scale) {
				float rows[2] = { i, i + scale };
				float cols[2] = { j, j + scale };

				count += check_bin(rows, cols, a, m);
			}
		}
		return count;
	}

	void linearRegression(float* coeffs, float* x_vals_inv, float* y_vals, int n_vals) {
		/* 
			simple linear regression taken from https://www.bragitoff.com/2015/09/c-program-to-linear-fit-the-data-using-least-squares-method/ 
			Note: proves to be untrustworthy in getting spectrum of fractal dimensions.
		*/
		float xsum = 0, x2sum = 0, ysum = 0, xysum = 0;

		for (int i = 0; i < n_vals; i++) {
			xsum = xsum + (x_vals_inv[i]);
			ysum = ysum + y_vals[i];
			x2sum = x2sum + pow((x_vals_inv[i]), 2);
			xysum = xysum + (x_vals_inv[i]) * y_vals[i];
		}

		coeffs[0] = (n_vals * xysum - xsum * ysum) / (n_vals * x2sum - xsum * xsum); /* slope */
		coeffs[1] = (x2sum * ysum - xsum * xysum) / (x2sum * n_vals - xsum * xsum); /* intercept */

	}

	float calcFracDim(int a, int m, int n_scales = 10, int n_offsets = 5, float min_exp = 2, float max_exp = 5) {
		/*
			Calculate fractal dimension values for a given array at
				number of scales n_scales
				number of offsets n_offests for count of boxes minimizin
				min_exp and max_exp define the scaling interval by box side length (these variables are exponents, not values)

				small-boxes scaling interval - 0.05 up to 2 (or 3)
				big-boxes scaling interval - 2 (or 3) up to quarter or half length of initial lattice side length
		
		*/
		float dim[2] = { 0 };
		float* scales = new float[n_scales]; /* input param */
		float* bin_vals = new float[n_scales];

		makeLogspace(scales, n_scales, min_exp, max_exp);
		

		for (int i = 0; i < n_scales; i++) {
			int* temp = new int[n_offsets];
			float* offsets = new float[n_offsets]; /* input param */

			makeLinspace(offsets, n_offsets, 0.01, scales[i] - 0.01);

			for (int j = 0; j < n_offsets; j++) {
				temp[j] = bin_count(scales[i], offsets[j], a, m);
			}

			bin_vals[i] = *min_element(temp, temp + n_offsets);

			delete[] offsets;
			delete[] temp;

			bin_vals[i] = log(bin_vals[i]);
			scales[i] = log((1 / scales[i]));
		}

		linearRegression(dim, scales, bin_vals, n_scales);

		delete[] bin_vals;
		delete[] scales;

		return dim[0];
	}

	/* 
		Functions below are necessary for mixed dynamics implementation
		Competing dynamics makes selection of dynamics and state (and spatial diversity) parameter calculation
		MC_Metropolis or MC_Kawasaki makes single dynamic implementation employing "sparse" matrix entries vectors:
			works very good when equilibriating Kawasaki type of dynamics (for small p)
			reduces speed of equilibriation at high p values. Metropolis step is very slow
		MC_Metropolis_clean or MC_Kawasaki_clean are single dynamics implementation without using vectors for mapping values

		Kawasaki dynamcis below are global
		Exponents are calculated for single dynamic calculation as they dont change values for one instance of object

		Equilibriation values need to be higher than usual (10^5 or 10^6)
		Averaging can be minimal (10 - 10^2) when lattice is equilibriated

		95% confidence gap may be calculated by sorting the Fdim or I averagin values arrays
	*/

	void competing_dynamics(float* ptr_E, float* ptr_M, float* ptr_M_up, float* ptr_M_down, float* ptr_I, float* ptr_I_up, float* ptr_I_down, float* ptr_F, float* ptr_F_up, float* ptr_F_down) { /* random number values will be passed from this function to MC step functions*/

		random_device rd; // obtain a random number from hardware
		mt19937 gen(rd()); // seed the generator		
		uniform_real_distribution<> distr1(0, 1); // define the range for rand [0-1)

		float exponents1[4] = { exp(-4 * beta1), exp(-8 * beta1), exp(-12 * beta1), exp(-16 * beta1) };
		float exponents2[4] = { exp(-4 * beta2), exp(-8 * beta2), exp(-12 * beta2), exp(-16 * beta2) };

		/* Equilibriation is for for same amount of iterations as averaging */
		for (int i = 0; i < t1; i++) {

			float probability = distr1(gen);

			if (probability <= dynamics_parameter) {
				MC_Metropolis_clean(exponents1, gen);
			}
			else {
				MC_Kawasaki_clean(exponents2, gen);
			}
		}
		printf("\t\t Start of averaging\n");
		for (int i = 0; i < t2; i++) {

			float probability = distr1(gen);

			if (probability <= dynamics_parameter ) { /* Let's say first is Metropolis */

				/* Calling Metropolis dynamics*/
				MC_Metropolis_clean(exponents1, gen);
			}
			else { /* Second is Kawasaki */

				/* Calling Kawasaki dynamics. They will not be executed if m is at either of limit values*/
				MC_Kawasaki_clean(exponents2, gen);
			}

			Mag[i] = calcMag() / (powf(N, 2));

			if (Mag[i] >= 0) {
				FracDim[i] = calcFracDim(1, Mag[i]);
			}
			else {
				FracDim[i] = calcFracDim(-1, Mag[i]);
			}

			Energy[i] = calcEnergy();
			Mag[i] = abs(Mag[i]); // Padariau kad butu absolute value
			//Index[i] = calcIndex();

		}
		
		for (int i = 1; i < t2; i++) {
			Energy[0] += Energy[i];
			Mag[0] += Mag[i];
			//Index[0] += Index[i];
			FracDim[0] += FracDim[i];
		}
		//sort(FracDim, FracDim + t2);
		*ptr_F = FracDim[0] / 2;
		//*ptr_F_down = FracDim[t2 / 20];
		//*ptr_F_up = FracDim[19 * t2 / 20];
		*ptr_E = Energy[0]/t2;
		*ptr_M = Mag[0]/t2;
		*ptr_I = Index[0]/t2;
		

		
	}

	void MC_Metropolis_clean(float* exponents, mt19937& gen) { /* need to implement spin swap*/

		uniform_int_distribution<> distrN(0, N - 1); // define the range for pos
		uniform_real_distribution<> distr(0, 1); // define the range for rand 0-1

		int n_steps = powf(N, 2);
		for (int i = 0; i < n_steps; i++) {
			int select_pos[2] = { distrN(gen), distrN(gen) };

			int dE = 2 * array[select_pos[0]][select_pos[1]] * stateNeighbours(select_pos);

			if (dE <= 0) {
				array[select_pos[0]][select_pos[1]] *= -1;

				/* transfer value in spin lists: */
				//transferSpin(select_pos); /* an already updated value is passed ! */
			}
			else if (distr(gen) < exponents[dE / 4 - 1]) {
				array[select_pos[0]][select_pos[1]] *= -1;
				//transferSpin(select_pos);
			}
		}
	}

	void MC_Kawasaki_clean(float* exponents, mt19937& gen) {

		uniform_int_distribution<> distr_dim(0, N - 1);
		uniform_real_distribution<> distr(0, 1); // define the range for rand [0-1)

		int n_steps = powf(N, 2);
		for (int i = 0; i < n_steps; i++) {

			int select_pos_up[2] = { distr_dim(gen), distr_dim(gen) };
			int select_pos_down[2] = { distr_dim(gen), distr_dim(gen) };

			if (array[select_pos_up[0]][select_pos_up[1]] != array[select_pos_down[0]][select_pos_down[1]]) {
				int dEu = array[select_pos_up[0]][select_pos_up[1]] * stateNeighbours(select_pos_up) + array[select_pos_down[0]][select_pos_down[1]] * stateNeighbours(select_pos_down);

				array[select_pos_up[0]][select_pos_up[1]] *= -1;
				array[select_pos_down[0]][select_pos_down[1]] *= -1;

				int dEv = array[select_pos_up[0]][select_pos_up[1]] * stateNeighbours(select_pos_up) + array[select_pos_down[0]][select_pos_down[1]] * stateNeighbours(select_pos_down);

				
				if (((dEu - dEv) <= 0) || (distr(gen) < exponents[(dEu - dEv) / 4 - 1])) {
				//if ((dEu - dEv) >= 0){
					//exchangeSpins(pos[0], pos[1]);
					int a;
				}
				else {
					array[select_pos_up[0]][select_pos_up[1]] *= -1;
					array[select_pos_down[0]][select_pos_down[1]] *= -1;
				}
			}

		}
	}

	void MC_Metropolis(float *exponents, mt19937& gen) { /* need to implement spin swap*/
		
		uniform_int_distribution<> distrN(0, N - 1); // define the range for pos
		uniform_real_distribution<> distr(0, 1); // define the range for rand 0-1

		int n_steps = powf(N, 2);
		for (int i = 0; i < n_steps; i++) {
			int select_pos[2] = { distrN(gen), distrN(gen) };

			int dE = 2 * array[select_pos[0]][select_pos[1]] * stateNeighbours(select_pos);

			if (dE <= 0) {
				array[select_pos[0]][select_pos[1]] *= -1;

				/* transfer value in spin lists: */
				transferSpin(select_pos); /* an already updated value is passed ! */
			}
			else if (distr(gen) < exponents[dE / 4 - 1]) {
				array[select_pos[0]][select_pos[1]] *= -1;
				transferSpin(select_pos);
			}
		}
	}

	void MC_Kawasaki(float *exponents, mt19937& gen) {

		if ((up_x.size() == 0) || (down_x.size() == 0)) {
			//cout << " One of the vectors is empty and therefore no Kawasaki dynamic is possible\n";
			return;
		}

		uniform_int_distribution<> distr_up(0, up_x.size() - 1); // two ranges for up and down vectors
		uniform_int_distribution<> distr_down(0, down_x.size() - 1); // two ranges for up and down vectors
		uniform_real_distribution<> distr(0, 1); // define the range for rand [0-1)

		int n_steps = powf(N, 2);
		for (int i = 0; i < n_steps; i++) {
			int pos[2] = { distr_up(gen), distr_down(gen) };
			int select_pos_up[2] = { up_x[pos[0]], up_y[pos[0]] };
			int select_pos_down[2] = { down_x[pos[1]], down_y[pos[1]] };

			int dEu = array[select_pos_up[0]][select_pos_up[1]] * stateNeighbours(select_pos_up) + array[select_pos_down[0]][select_pos_down[1]] * stateNeighbours(select_pos_down);

			array[select_pos_up[0]][select_pos_up[1]] *= -1;
			array[select_pos_down[0]][select_pos_down[1]] *= -1;

			int dEv = array[select_pos_up[0]][select_pos_up[1]] * stateNeighbours(select_pos_up) + array[select_pos_down[0]][select_pos_down[1]] * stateNeighbours(select_pos_down);

			if (((dEu - dEv) <= 0) || (distr(gen) < exponents[(dEu - dEv) / 4 - 1])) {
				exchangeSpins(pos[0], pos[1]);
			}
			else {
				array[select_pos_up[0]][select_pos_up[1]] *= -1;
				array[select_pos_down[0]][select_pos_down[1]] *= -1;
			}

		}
	}
};

void printResults(dataIsing results, string filenames);

void makeLinspaceGlobal(float* arr);

void printEquilibrium(float* E, float* M, float* I, int* P, float* F, int num);

int main() {

	/*
		The code below is for equilibriation implementation
		It was only used for initial tests of code
	*/
	/*
	compEquilib data;
	data.N = 128;
	data.T1 = 1.5;
	data.T2 = 1.5;
	data.p = 0.1;
	data.m = 0.5;
	data.uni = false;

	for (int i = 0; i < 9; i++) {
		cout <<"p = 0."<< i+1 << " started" << endl;
		float* E = new float[iterations];
		float* M = new float[iterations];
		float* I = new float[iterations];
		int* P = new int[iterations];
		float* F = new float[iterations];

		data.p = data.p + 0.1 * i;
		Competing model(data.N, data.T1, data.T2, data.p, data.m, iterations, true, true, data.uni);
		model.competing_dynamics(E, M, I, P);
		printEquilibrium(E, M, I, P, F, i+1);

		delete[] E;
		delete[] M;
		delete[] I;
		delete[] P;
		delete[] F;
	}
	*/
	
	/*
		Below is the main function for dynamics over a range of parameters T1, T2 and p calculations
		We set number of p values and initiate in the for cycle
		T1 (Metropolis temperature) values are set by constants in the beggining of .cpp file
		T2 values should be set inside the cycle for various cases of Kawasaki dynamics
	*/
	const int num_points = 4;
	float p[num_points] = { 0 };
	float m = 0;
	int it2 = 100;
	int it1 = 100000;
	
	dataIsing results;


	for (int j = 0; j < num_points; j++) {
		p[j] = 0.4 + 0.2 * j;
		printf("Start of p=%0.2f\n", p[j]);

		dataIsing results;
		//makeLinspaceGlobal(results.T1);

		ifstream file;
		const string fin = "T-input.txt";
		file.open(fin);

		for (int i = 0; i < data_points; i++) {
			file >> results.T1[i];
		}

		file.close();

		for (int k = 0; k < data_points; k++) {
			results.T2[k] = 0;
		}

		string filenames = filepath + to_string(int(10*p[j])) + ending;
		//cout << filenames << endl;

		
		for (int i = 0; i < data_points; i++) {
			printf("\tStart of data point calculations for number %10d \n", i);
			Competing model(64, results.T1[i], results.T2[i], p[j], m, it1, it2, false);
			model.competing_dynamics(&results.E[i], &results.M[i], &results.M_up[i], &results.M_down[i], &results.I[i], &results.I_95up[i], &results.I_95down[i], &results.F[i], &results.F_95up[i], &results.F_95down[i]);
		}
		
		printResults(results, filenames);
	}
	return 0;
}

void printEquilibrium(float *E, float *M, float *I, int *P, float *F, int num) {

	ofstream file;
	file.open(filename+to_string(num)+".txt");

	for (int i = 0; i < iterations; i++) {
		file << E[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < iterations; i++) {
		file << M[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < iterations; i++) {
		file << I[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < iterations; i++) {
		file << P[i] << "\t";
	}

	for (int i = 0; i < iterations; i++) {
		file << F[i] << "\t";
	}
	file << "\n";

	file.close();
}

void printResults(dataIsing results, string filenames) {

	ofstream file;
	file.open(filenames);

	for (int i = 0; i < data_points; i++) {
		file << results.T1[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.T2[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.E[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.M[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.M_up[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.M_down[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.I[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.I_95up[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.I_95down[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.F[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.F_95up[i] << "\t";
	}
	file << "\n";

	for (int i = 0; i < data_points; i++) {
		file << results.F_95down[i] << "\t";
	}
	file << "\n";

	file.close();
}

void makeLinspaceGlobal(float* arr) {
	float diff = (max_val - min_val) / (data_points - 1);
	for (int i = 0; i < data_points; i++) {
		arr[i] = min_val + diff * i;
	}
}