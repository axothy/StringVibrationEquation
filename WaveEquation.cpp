#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <array>
#include <stdint.h>

#define _USE_MATH_DEFINES

#if defined(_USE_MATH_DEFINES) && !defined(_MATH_DEFINES_DEFINED)
#define _MATH_DEFINES_DEFINED
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076

#endif  /* _USE_MATH_DEFINES */


#define N 11 // 0 < x < 1
#define K 51 //0 < t < 0.5

enum ListOfMediums {
	DryAir = 0,
	Hydrogen = 1,
	Nitrogen = 2,
	Helium = 3,
	Oxygen = 4,
	CarbonDioxide = 5,

	Acetone = 6,
	Water = 7,
	Petrol = 8,
	Vodka = 9,

	One = 10 //by exercise c = 1
};

struct MediumData {
	const int Medium_ID;
	const double c_SoundWaveSpeed;
};

namespace MediumDataBase {
	std::vector<MediumData> mediums = {
		{0, 331.0},
		{1, 1284.0},
		{2, 334.0},
		{3, 955.0},
		{4, 316.0},
		{5, 259.0},

		{6, 1192.0},
		{7, 1460.0},
		{8, 1170.0},
		{9, 1180.0},

		{10, 1}
	};
}

struct IntegrateParameters {
	double delta_tStep;
	double h_xStep;
};

class StringVibrationSolver {

private:
	inline double InitialValueFunction(double x);
	inline double DerivativeIVF(double x);
	inline double LeftBoundaryConditionFunction(double t);
	inline double RightBoundaryConditionFunction(double t);

public:
	void SetCoordinateXGrid(std::array<double, N>& x, const IntegrateParameters& params);
	void SetTimeGrid(std::array<double, K>& t, const IntegrateParameters& params);

	void getSolutionExplicit(const MediumData& medium, double(*U)[N], std::array<double, N>& x, 
		std::array<double, K>& t, const IntegrateParameters& params);
	void getSolutionImplicit(const MediumData& medium, double(*U)[N], std::array<double, N>& x, 
		std::array<double, K>& t, const IntegrateParameters& params);
};


//НУ T(x,0) 
inline double StringVibrationSolver::InitialValueFunction(double x)
{
	return (x * x + 0.5) * cos(M_PI * x);
}

//derived NU
inline double StringVibrationSolver::DerivativeIVF(double x)
{
	return (x + 0.7) * (x + 0.7);

}

//ГУ1 T(0, t)
inline double StringVibrationSolver::LeftBoundaryConditionFunction(double t)
{
	return 0.5;
}

//ГУ2 T(l,t)
inline double StringVibrationSolver::RightBoundaryConditionFunction(double t)
{
	return 2 * t - 1.5;
}

//Создание сетки по координате х
void StringVibrationSolver::SetCoordinateXGrid(std::array<double, N>& x, const IntegrateParameters& params)
{
	x[0] = 0;
	for (uint8_t i = 1; i < N; ++i) {
		x[i] = i * params.h_xStep;
	}
}

void StringVibrationSolver::SetTimeGrid(std::array<double, K>& t, const IntegrateParameters& params)
{
	t[0] = 0;
	for (uint8_t k = 1; k < K; ++k) {
		t[k] = k * params.delta_tStep;
	}
}

void StringVibrationSolver::getSolutionExplicit(const MediumData& medium, double(*U)[N], std::array<double, N>& x,
	std::array<double, K>& t, const IntegrateParameters& params)
{

	for (uint8_t i = 0; i < N; ++i)
	{
		U[0][i] = InitialValueFunction(x[i]);
		U[1][i] = InitialValueFunction(x[i]) + DerivativeIVF(x[i]) * params.delta_tStep;
	}


	for (uint8_t k = 0; k < K; ++k)
	{
		U[k][0] = LeftBoundaryConditionFunction(t[k]);
		U[k][N-1] = RightBoundaryConditionFunction(t[k]);
	}

	/*Для сокращения выражений*/
	double h = params.h_xStep;
	double dt = params.delta_tStep;
	double c = medium.c_SoundWaveSpeed;

	/*Решение по явной схеме*/
	for (uint8_t k = 1; k < K-1; ++k)
	{
		for (uint8_t i = 1; i < N - 1; ++i)
		{
			U[k + 1][i] = ((c * c) * (dt * dt) / (h * h)) * (U[k][i + 1] - 2 * U[k][i] + U[k][i - 1]) +
				2 * U[k][i] - U[k - 1][i];
		}
	}
}

void StringVibrationSolver::getSolutionImplicit(const MediumData& medium, double(*U)[N], std::array<double, N>& x, std::array<double, K>& t,
	const IntegrateParameters& params)
{
	for (uint8_t i = 0; i < N; ++i)
	{
		U[0][i] = InitialValueFunction(x[i]);
		U[1][i] = InitialValueFunction(x[i]) + DerivativeIVF(x[i]) * params.delta_tStep;
	}


	for (uint8_t k = 0; k < K; ++k)
	{
		U[k][0] = LeftBoundaryConditionFunction(t[k]);
		U[k][N - 1] = RightBoundaryConditionFunction(t[k]);
	}


	double h = params.h_xStep;
	double dt = params.delta_tStep;
	double c = medium.c_SoundWaveSpeed;
	
	double A = 1.0 / (h * h);
	double C = 1.0 / (h * h);
	double B = ((h*h) + 2 * (c*c) * (dt*dt)) / ((h*h) * (c*c) * (dt*dt));
	std::array<double, N> F, P, Q;

	for (uint8_t k = 1; k < K - 1; ++k)
	{

		for (uint8_t i = 0; i < N; i++)
		{
			F[i] = -((U[k - 1][i]) / (c * c * dt * dt)) + U[k][i] * (2 / (c * c * dt * dt));
		}
		P[0] = C / B;
		Q[0] = F[0] / B;

		for (uint8_t j = 1; j < N; ++j)
		{
			P[j] = C / (B - A * P[j - 1]);
			Q[j] = (F[j] + A * Q[j - 1]) / (B - A * P[j - 1]);
		}

		for (uint8_t j = N - 2; j > 0; --j)
		{
			U[k + 1][j] = P[j] * U[k + 1][j + 1] + Q[j];
		}
	}
	

}


int main()
{
	using namespace MediumDataBase;

	StringVibrationSolver mySolver;
	IntegrateParameters params;

	params.h_xStep = 0.1;
	params.delta_tStep = 0.01;

	double U[K][N] = { 0 }; 
	std::array<double, N> x;
	std::array<double, K> t;

	mySolver.SetCoordinateXGrid(x, params);
	mySolver.SetTimeGrid(t, params);

	/*Здесь можно задать материал (см. enum)*/
	mySolver.getSolutionImplicit(mediums[One], U, x, t, params);

	for (int i = 0; i < N; i++)
		printf("%lf ", x[i]);
	printf("\n\n");
	for (int k = 0; k < K; k++)
		printf("%lf ", t[k]);
	printf("\n\n");

	//$t_{ 0 } = 0$ & 0.900000 & 1.080000 & 1.220000 & 1.320000 & 1.380000 & 1.400000 & 1.380000 & 1 & 1 & 1 & 1 \\ \hline

	for (int i = 0; i < K; ++i)
	{
		printf("$t_{ %d } = %.2lf$ ", i, t[i]);
		for (int j = 0; j < N; ++j)
			printf("& %lf\t", U[i][j]);
		printf("\\\\  \\hline");
		printf("\n\n");
	}

	return 0;
}