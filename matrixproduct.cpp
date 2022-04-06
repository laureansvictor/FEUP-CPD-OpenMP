#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <papi.h>
#include <math.h>
#include <omp.h>

using namespace std;

#define SYSTEMTIME clock_t

int EventSet;

double OnMult(int m_ar, int m_br)
{

	SYSTEMTIME Time1, Time2;

	char st[100];
	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);

	// Start counting
	int ret = PAPI_start(EventSet);
	if (ret != PAPI_OK)
		cout << "ERROR: Start PAPI" << endl;

	Time1 = clock();
	#pragma omp parallel for
	for (i = 0; i < m_ar; i++)
	{
		for (j = 0; j < m_br; j++)
		{
			temp = 0.0;
			for (k = 0; k < m_ar; k++)
			{
				temp += pha[i * m_ar + k] * phb[k * m_br + j];
			}
			phc[i * m_ar + j] = temp;
		}
	}

	Time2 = clock();
	sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
	cout << st;

	// display 10 elements of the result matrix to verify correctness
	cout << "Result matrix: " << endl;
	for (i = 0; i < 1; i++)
	{
		for (j = 0; j < min(10, m_br); j++)
			cout << phc[j] << " | ";
	}
	cout << endl;

	free(pha);
	free(phb);
	free(phc);

	return (double)(Time2 - Time1) / CLOCKS_PER_SEC;
}

double OnMultLine(int m_ar, int m_br)
{
	SYSTEMTIME Time1, Time2;

	char st[100];
	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phc[i * m_br + j] = (double)(0);

	// Start counting
	int ret = PAPI_start(EventSet);
	if (ret != PAPI_OK)
		cout << "ERROR: Start PAPI" << endl;
	Time1 = clock();
    
	#pragma omp parallel
	for (i = 0; i < m_ar; i++)
	{
		for (k = 0; k < m_ar; k++)
		{
			for (j = 0; j < m_br; j++)
			{
				phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];
			}
		}
	}

	Time2 = clock();
	sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
	cout << st;

	cout << "Result matrix: " << endl;
	for (i = 0; i < 1; i++)
	{
		for (j = 0; j < min(10, m_br); j++)
			cout << phc[j] << " | ";
	}
	cout << endl;

	free(pha);
	free(phb);
	free(phc);

	return (double)(Time2 - Time1) / CLOCKS_PER_SEC;
}

double OnMultBlock(int m_ar, int m_br, int blk_sizes)
{
	SYSTEMTIME Time1, Time2;

	char st[100];
	double temp;
	int i, j, k;

	double *pha, *phb, *phc;

	pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
	phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

	for (i = 0; i < m_ar; i++)
		for (j = 0; j < m_ar; j++)
			pha[i * m_ar + j] = (double)1.0;

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phb[i * m_br + j] = (double)(i + 1);

	for (i = 0; i < m_br; i++)
		for (j = 0; j < m_br; j++)
			phc[i * m_br + j] = (double)(0);
	// Start counting
	int ret = PAPI_start(EventSet);
	if (ret != PAPI_OK)
		cout << "ERROR: Start PAPI" << endl;
	Time1 = clock();

	for (i = 0; i < m_ar; i += blk_sizes)
	{
		for (j = 0; j < m_br; j += blk_sizes)
		{
			for (int blk_row = i; blk_row < min(i + blk_sizes, m_ar); blk_row++)
			{
				for (k = 0; k < m_ar; k++)
				{
					for (int blk_col = j; blk_col < min(j + blk_sizes, m_ar); blk_col++)
					{
						phc[blk_row * m_ar + blk_col] += pha[blk_row * m_ar + k] * phb[k * m_br + blk_col];
					}
				}
			}
		}
	}

	Time2 = clock();
	sprintf(st, "Time: %3.3f seconds\n", (double)(Time2 - Time1) / CLOCKS_PER_SEC);
	cout << st;

	cout << "Result matrix: " << endl;
	for (i = 0; i < 1; i++)
	{
		for (j = 0; j < min(10, m_br); j++)
			cout << phc[j] << " | ";
	}
	cout << endl;

	free(pha);
	free(phb);
	free(phc);

	return (double)(Time2 - Time1) / CLOCKS_PER_SEC;
}

float produtoInterno(float *v1, float *v2, int col)
{
	int i;
	float soma = 0.0;

	for (i = 0; i < col; i++)
		soma += v1[i] * v2[i];

	return (soma);
}

void handle_error(int retval)
{
	printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
	exit(1);
}

void init_papi()
{
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT && retval < 0)
	{
		printf("PAPI library version mismatch!\n");
		exit(1);
	}
	if (retval < 0)
		handle_error(retval);

	std::cout << "PAPI Version Number: MAJOR: " << PAPI_VERSION_MAJOR(retval)
			  << " MINOR: " << PAPI_VERSION_MINOR(retval)
			  << " REVISION: " << PAPI_VERSION_REVISION(retval) << "\n";
}

// int setupPAPI(int &EventSet)
// {
// 	int _events[NUMEVENTS] = {PAPI_DP_OPS, PAPI_L2_DCM, PAPI_L2_DCA, PAPI_L2_TCM, PAPI_L2_TCA, PAPI_TOT_CYC};
// 	int ret = PAPI_library_init(PAPI_VER_CURRENT);
// 	if (ret != PAPI_VER_CURRENT)
// 		std::cout << "FAIL" << endl;

// 	ret = PAPI_create_eventset(&EventSet);
// 	if (ret != PAPI_OK)
// 		cout << "ERROR: create eventset" << endl;

// 	ret = PAPI_add_events(EventSet, _events, 6);
// 	if (ret != PAPI_OK)
// 		cout << "ERROR: Starting Papi Events" << endl;
// 	return ret;
// }

void test(int line, int end, int step, fstream &logfile, double (*func)(int, int))
{
	double time;
	int ret;
	long long values[4];

	logfile << "Line Multiplication: from " << line << " to " << end << " step:" << step << endl;
	logfile << "lineCol"
			<< " | "
			<< "time"
			<< " | "
			<< "CYCLES"
			<< " | "
			<< "L2 Miss ratio"
			<< " | "
			<< "L2 Hit ratio"
			<< " | "
			<< "GFLOPS"
			<< " | "
			<< "Inst per Cycle"
			<< " | " << endl;

	for (int lineCol = line; lineCol <= end; lineCol += step)
	{
		time = func(lineCol, lineCol);

		ret = PAPI_stop(EventSet, values);

		printf("retorno da função de multiplicação %d\n", ret);
		if (ret != PAPI_OK)
			cout << "ERROR: Stop PAPI" << endl;
		double l2dmiss = (double)values[1] / ((double)values[1] + (double)values[0]);
		double gflops = (double)(2 * (pow(lineCol * lineCol, 3)) / time);
		double ipc = (double)values[3] / (double)values[2];
		printf("CYCLES: %lld \n", values[2]);
		printf("L2 DCA: %lld \n", values[0]);
		printf("L2 DCM: %lld \n", values[1]);
		printf("Total Instructions: %lld \n", values[3]);
		printf("%lld L2 Data Cache Misses (%.4lf%% misses) with %lld L2 Data Cache Access in %lld Cycles\n",
			   values[1], l2dmiss, values[0], values[2]);
		printf("GFLOPS %.10lf 10e9\n", gflops * 1 * 10e-9);
		printf("Inst per cycle %.10lf \n", ipc);

		ret = PAPI_reset(EventSet);
		if (ret != PAPI_OK)
			std::cout << "FAIL reset" << endl;

		logfile << lineCol << " | " << time << " | " << values[2] << " | " << l2dmiss << " | " << 1 - l2dmiss << " | " << gflops * 1 * 10e-9 << " | " << ipc << " | " << endl;
	}
}

void blockTest(int line, int end, int step, int blockSize, fstream &logfile)
{
	double time;
	int ret;
	// int EventSet = PAPI_NULL;
	long long values[4];
	// int ret = setupPAPI(EventSet);

	logfile << "Block Multiplication: from " << line << " to " << end << " step:" << step << endl;
	logfile << "lineCol"
			<< " | "
			<< "BlockSize"
			<< " | "
			<< "time"
			<< " | "
			<< "CYCLES"
			<< " | "
			<< "L2 Miss ratio"
			<< " | "
			<< "L2 Hit ratio"
			<< " | "
			<< "GFLOPS"
			<< " | "
			<< "Inst per Cycle"
			<< " | " << endl;

	for (int lineCol = line; lineCol <= end; lineCol += step)
	{
		time = OnMultBlock(lineCol, lineCol, blockSize);
		ret = PAPI_stop(EventSet, values);
		printf("retorno da função de multiplicação %d\n", ret);
		if (ret != PAPI_OK)
			cout << "ERROR: Stop PAPI" << endl;
		double l2dmiss = (double)values[1] / ((double)values[1] + (double)values[0]);
		double gflops = (double)(2 * (pow(lineCol * lineCol, 3)) / time);
		double ipc = (double)values[3] / (double)values[2];
		printf("CYCLES: %lld \n", values[2]);
		printf("L2 DCA: %lld \n", values[0]);
		printf("L2 DCM: %lld \n", values[1]);
		printf("Total Instructions: %lld \n", values[3]);
		printf("%lld L2 Data Cache Misses (%.4lf%% misses) with %lld L2 Data Cache Access in %lld Cycles\n",
			   values[1], l2dmiss, values[0], values[2]);
		printf("GFLOPS %.10lf 10e9\n", gflops * 1 * 10e-9);
		printf("Inst per cycle \n%.10lf", ipc);

		ret = PAPI_reset(EventSet);
		if (ret != PAPI_OK)
			std::cout << "FAIL reset" << endl;

		logfile << lineCol << " | " << blockSize << " | " << time << " | " << values[2] << " | " << l2dmiss << " | " << 1 - l2dmiss << " | " << gflops * 1 * 10e-9 << " | " << ipc << " | " << endl;
	}
}

int main(int argc, char *argv[])
{

	char c;
	int lin, col, nt = 1;
	int op;
	int blockSize;
	double time = 0;

	int EventSet = PAPI_NULL;
	long long values[4];
	int ret;

	fstream logfile;
	logfile.open("logfile.txt", ifstream::out);

	int _events[4] = {PAPI_L2_DCA, PAPI_L2_DCM, PAPI_TOT_CYC, PAPI_TOT_INS};

	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
	{
		fprintf(stderr, "PAPI Library initialization error! %d\n", __LINE__);
		exit(1);
	}

	ret = PAPI_create_eventset(&EventSet);
	if (ret != PAPI_OK)
		cout << "ERROR: create eventset" << endl;
	ret = PAPI_add_events(EventSet, _events, 4);
	if (ret != PAPI_OK)
		printf("ERROR ADDING EVENT: A total of %d events were added", ret);

	op = 1;
	do
	{
		cout << endl
			 << "1. Multiplication" << endl;
		cout << "2. Line Multiplication" << endl;
		cout << "3. Block" << endl;
		cout << "4. Multiplication test 600 -> 3000, step: 400" << endl;
		cout << "5. Line Multiplication test 600 -> 3000, step: 400" << endl;
		cout << "6. Line Multiplication test 4096 -> 10240, step: 2048" << endl;
		cout << "7. Block Multiplication test 4096 -> 10240, step: 2048, with blocks 128, 256 and 512" << endl;
		cout << "9. All rests" << endl;
		cout << "Selection?: ";
		cin >> op;
		if (op == 0)
			break;
		if (op == 1 || op == 2 || op == 3)
		{
			printf("Dimensions: lins=cols ? ");
			cin >> lin;
			col = lin;
		}

		switch (op)
		{
		case 1:
			time = OnMult(lin, col);
			break;
		case 2:
			time = OnMultLine(lin, col);
			break;
		case 3:
			cout << "Block Size?";
			cin >> blockSize;
			time = OnMultBlock(lin, col, blockSize);
			break;
		case 4:
			test(600, 3000, 400, logfile, &OnMult);
			break;
		case 5:
			test(600, 3000, 400, logfile, &OnMultLine);
			break;
		case 6:
			test(2048, 10240, 2048, logfile, &OnMultLine);
			break;
		case 7:
			blockTest(4096, 10240, 2048, 128, logfile);
			blockTest(4096, 10240, 2048, 256, logfile);
			blockTest(4096, 10240, 2048, 512, logfile);
			break;
		case 9:
			test(600, 3000, 400, logfile, &OnMult);
			test(600, 3000, 400, logfile, &OnMultLine);
			test(2048, 10240, 2048, logfile, &OnMultLine);
			blockTest(4096, 10240, 2048, 128, logfile);
			blockTest(4096, 10240, 2048, 256, logfile);
			blockTest(4096, 10240, 2048, 512, logfile);
			break;
		}

		ret = PAPI_stop(EventSet, values);
		printf("retorno da função de multiplicação %d\n", ret);
		if (ret != PAPI_OK)
			cout << "ERROR: Stop PAPI" << endl;

		double l2dmiss = (double)values[1] / ((double)values[1] + (double)values[0]);
		double gflops = (double)(2 * (pow(lin * col, 3)) / time);
		double ipc = (double)values[3] / (double)values[2];
		printf("CYCLES: %lld \n", values[2]);
		printf("L2 DCA: %lld \n", values[0]);
		printf("L2 DCM: %lld \n", values[1]);
		printf("Total Instructions: %lld \n", values[3]);
		printf("%lld L2 Data Cache Misses (%.4lf%% misses) with %lld L2 Data Cache Access in %lld Cycles\n",
			   values[1], l2dmiss, values[0], values[2]);
		printf("GFLOPS %.10lf \n", gflops * 1 * 10e-9);
		printf("Inst per cycle %.10lf \n", ipc);

		ret = PAPI_reset(EventSet);
		if (ret != PAPI_OK)
			std::cout << "FAIL reset" << endl;

	} while (op != 0);

	ret = PAPI_remove_events(EventSet, _events, 4);
	if (ret != PAPI_OK)
		std::cout << "FAIL remove events" << endl;

	ret = PAPI_destroy_eventset(&EventSet);
	if (ret != PAPI_OK)
		std::cout << "FAIL destroy" << endl;
}
