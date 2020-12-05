//#define Debug

#include <iostream>
#include <complex>
#include <time.h>
#include <vector>
#include <thread>
#include <cstdlib>
#include <omp.h>

using namespace::std;
using namespace::chrono;

//unsigned int size_ = 1024;

unsigned int size_ = 1024;
unsigned int sizeSquared = size_*size_;

// 1 замер - > 
// 2 замер - > 
// 3 замер - >
//
//
//
//

//Вывод матрицы
void PrintMatrix(complex<double>* arr)
{
	for (unsigned int i = 0; i < size_; i++)
	{
		for (unsigned int j = 0; j < size_; j++)
		{
			cout << arr[i * size_ + j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
}

//Возврат рандомного значения типа double в указанном диапазоне
double randomNumber(int min, int max)
{
	//1.0 + 4.0 * rand() / (double)RAND_MAX;
	return (double)rand() / (RAND_MAX + 1) * (double)(max - min) + min;
	//return min + rand() % (1000 * (max - min)) / 1000.0;
}

//Оценка сложности алгоритма
unsigned long long algorithmСomplexity(unsigned int n)
{
	return (unsigned long long)(2 * pow(n, 3));
}

//Оценка производительности
double performance(double t, unsigned long long c)
{
	return (((double)c / t) * (double)pow(10, -6));
}

complex<double>* transpose(complex<double>* matrix) 
{
	complex<double> t;
	for (int i = 0; i < size_; ++i)
	{
		for (int j = i; j < size_; ++j)
		{
			t = matrix[i*size_+j];
			matrix[i*size_+j] = matrix[j*size_+i];
			matrix[j*size_+i] = t;
		}
	}
	return matrix;
}

//========================================  ПЕРВЫЙ МЕТОД ПРОИЗВЕДЕНИЯ  ===========================================

complex<double>* firstMethod(const complex<double>* A, complex<double>* B)//Первый метод||по формуле из линейной алгебры
{
	complex <double>* Arr = new complex<double>[sizeSquared];
	complex <double>* tB = transpose(B);

	omp_set_num_threads(2000);
	{
		unsigned int i, j, k;
#pragma omp parallel for;
		{
			for (i = 0; i < size_; i++)
			{
				for (j = 0; j < size_; j++)
				{
					complex<double> temp = Arr[i * size_ + j];
					for (k = 0; k < size_; k++)
					{
						temp += A[i * size_ + k] * tB[j * size_ + k];
					}
					Arr[i * size_ + j] = temp;
				}
			}
		}
		/*for(i = 0; i < size_; i++)
		{
			for (j = 0; j < size_; j++)
			{
				complex<double>tmp = (0, 0);
				for (k = 0; k < size_; k++)
				{
					tmp += A[i * size_ + k] * tB[j * size_ + k];
				}
				Arr[i * size_ + j] = tmp;
			}
		}*/
		/*for ( int i = 0; i < size_; i++)
		{
			for ( int j = 0; j < size_; j++)
			{
				for ( int k = 0; k < size_; k++)
				{
					Arr[i][j] += A[i][k] * B[k][j];
				}
			}
		}*/
	}
	
	return Arr;
}
void firstAlgorithm(const complex<double>* Array1, complex<double>* Array2)//Первый алгоритм
{
	cout << "Умножаем матрицы 1 способом: " << endl;

	auto startTime1 = high_resolution_clock::now();//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	complex <double>* Arr1 = firstMethod(Array1, Array2);

	auto endTime1 = high_resolution_clock::now();//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	duration<float> Time1 = endTime1 - startTime1;

	unsigned long long Complexity1 = algorithmСomplexity(size_);//Сложность 1 алгоритма
	cout << "Время произведения матриц первым методом: " << Time1.count() << "с" << endl;
	cout << "Производительность первого алгоритма произведения матриц = " << performance((double)Time1.count(), Complexity1) << " MFlops" << endl;

#if defined(Debug)//Вывод
	PrintMatrix(Arr1);
#endif // gebug
	system("pause");
}

//========================================  ВТОРОЙ МЕТОД ПРОИЗВЕДЕНИЯ  ===========================================

void secondMethod(vector<vector<complex<double>>> arr1, vector<vector<complex<double>>> arr2)//Второй метод||результат работы функции cblas_zgemm из библиотеки BLAS
{

}

void thirdMethod(vector<vector<complex<double>>> arr1, vector<vector<complex<double>>> arr2)//Третий метод||оптимизированный алгоритм по выбору
{

}

//============================================================================================= Конец  написания 3 алгоритмов умножения матриц


int main()
{
	setlocale(LC_ALL, "Russian");
	srand(time(0));


	cout << "Создаем двумерные массивы комплексных чисел! Ожидайте..." << endl;
	auto startTime = high_resolution_clock::now();//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Начальное время

	//Выделяем память под двумерные массивы
	complex<double>* Array1 = new complex<double>[sizeSquared];
	complex<double>* Array2 = new complex<double>[sizeSquared];

	//Инициализируем массивы комплексными числами
	for (unsigned int i = 0; i < size_; i++)
	{
		for (unsigned int j = 0; j < size_; j++)
		{
			/*Array1[i * size_ + j] = complex<double>(randomNumber(10, 100), randomNumber(10, 100));
			Array2[i * size_ + j] = complex<double>(randomNumber(10, 100), randomNumber(10, 100));*/
			Array1[i * size_ + j] = complex<double>(i,j);
			Array2[i * size_ + j] = complex<double>(j,i);
		}
	}
	auto endTime = high_resolution_clock::now();//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Конечное время
	duration<float> searchTime = endTime - startTime;

	unsigned long long c = algorithmСomplexity(size_);

	system("cls");
	cout << "Время создания и заполнения двух массивов: " << searchTime.count() << "с" << endl;
	cout << "Сложность алгоритма = " << c << endl;
	cout << "Производительность алгоритма = " << performance((double)searchTime.count(), c) << " MFlops" << endl << endl;


#if defined(Debug)//Вывод
	PrintMatrix(Array1);
	PrintMatrix(Array2);
#endif // gebug

	firstAlgorithm(Array1, Array2);//Из линейной алгебры

}


