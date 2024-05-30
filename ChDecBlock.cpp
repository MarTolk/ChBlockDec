#include<iostream>
#include<vector>
#include<cmath>


//========================================================================================================================
void PrMatrix(std::vector<std::vector<float>>& matrix)
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[i].size(); j++)
		{
			std::cout << matrix[i][j] << ' ';
		}
		std::cout << "\n";
	}
}
void PrVec(std::vector<float>& Vec)
{
	for (int i = 0; i < Vec.size(); i++)
	{
		std::cout << Vec[i] << ' ';
	}
	std::cout << "\n";
}
void WithZeroVec(std::vector<float>& Vec)
{
	for (int i = 0; i < Vec.size(); i++)
	{
		Vec[i] = 0.0;
	}
}
//========================================================================================================================
void SplitToBlock(std::vector<std::vector<float>>& A, std::vector<float>& Inp, int& MSide, int& BlSide)
{
	std::vector<float>Vec;
	Vec.resize(MSide * MSide);

	int VertBlNumber = 0, HorizBlockNumber = 0, ElInBlNumber = 0;
	int VertBlIndex = 0, HorizBlIndex = 0;
	int iBlock = 0, jBlock = 0;
	int BlockIndex = 0;
	int count = 0;

	VertBlNumber = MSide / BlSide;
	HorizBlockNumber = MSide / BlSide;
	ElInBlNumber = BlSide * BlSide;

	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			VertBlIndex = i / BlSide;
			HorizBlIndex = j / BlSide;
			iBlock = i % BlSide;
			jBlock = j % BlSide;
			BlockIndex = VertBlIndex * HorizBlockNumber + HorizBlIndex;
			Vec[(BlockIndex * ElInBlNumber) + (iBlock * BlSide + jBlock)] = Inp[i * MSide + j];
		}
	}
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			A[i][j] = Vec[count];
			count++;
		}
	}
	Vec.clear();
}

void CombBlock(std::vector<std::vector<float>>& Matrix, std::vector<float>& Vec, int& MSide, int& BlSide)
{
	std::vector<float>Vec1;
	Vec1.resize(MSide * MSide);
	int count = 0;
	for (int i = 0; i < Matrix.size(); i++)
	{
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			Vec1[count] = Matrix[i][j];
			count++;

		}
	}
	int VertBlNumber = 0, HorizBlockNumber = 0, ElInBlNumber = 0;
	VertBlNumber = MSide / BlSide;
	HorizBlockNumber = MSide / BlSide;
	ElInBlNumber = BlSide * BlSide;
	int VertBlIndex = 0, HorizBlIndex = 0;
	int iBlock = 0, jBlock = 0;
	int BlockIndex = 0;

	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			VertBlIndex = i / BlSide;
			HorizBlIndex = j / BlSide;
			iBlock = i % BlSide;
			jBlock = j % BlSide;
			BlockIndex = VertBlIndex * HorizBlockNumber + HorizBlIndex;
			Vec[i * MSide + j] = Vec1[(BlockIndex * ElInBlNumber) + (iBlock * BlSide + jBlock)];
		}
	}
	Vec1.clear();
}
void ChDecOneBl(std::vector<float>& SubA, std::vector<float>& SubL, int& BlSide)
{
	float Sum = 0.0;
	for (int i = 0; i < BlSide; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			Sum = 0.0;
			if (j == i)
			{
				for (int k = 0; k < j; k++)
				{
					Sum += std::pow(SubL[j * BlSide + k], 2);
				}
				SubL[j * BlSide + j] = std::sqrt(SubA[j * BlSide + j] - Sum);
			}
			else
			{
				for (int k = 0; k < j; k++)
				{
					Sum += (SubL[i * BlSide + k] * SubL[j * BlSide + k]);
				}
				SubL[i * BlSide + j] = (SubA[i * BlSide + j] - Sum) / SubL[j * BlSide + j];
			}
		}
	}
}
void CalcMat(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& L,
	std::vector<float>& MatOut, int& BlSide, int& IndI, int& IndJ, int& MTileSide)
{
	std::vector<float>Mult, Sum;
	Mult.resize(BlSide * BlSide);
	Sum.resize(BlSide * BlSide);

	for (int m = 0; m < IndJ; m++)
	{
		for (int i = 0; i < BlSide; i++)
		{
			for (int j = 0; j < BlSide; j++)
			{
				for (int k = 0; k < BlSide; k++)
				{
					Mult[i * BlSide + j] = Mult[i * BlSide + j] + (L[IndI * MTileSide + m][i * BlSide + k] * L[IndJ * MTileSide + m][j * BlSide + k]);
				}
				Sum[i * BlSide + j] = Sum[i * BlSide + j] + Mult[i * BlSide + j];
			}
		}
		WithZeroVec(Mult);
	}
	for (int i = 0; i < BlSide; i++)
	{
		for (int j = 0; j < BlSide; j++)
		{
			MatOut[i * BlSide + j] = A[IndI * MTileSide + IndJ][i * BlSide + j] - Sum[i * BlSide + j];
		}
	}
	Mult.clear();
	Sum.clear();
}
void LowInvTr(std::vector<float>& LowTr, std::vector<float>& Prom, std::vector<float>& MatOut, int& MSide)
{
	float PrElem = 0.0;
	std::vector<float>LInv;
	LInv.resize(MSide * MSide);
	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			if (i == j)
			{
				LInv[i * MSide + i] = 1 / LowTr[i * MSide + i];
			}
			else if (j < i)
			{
				for (int k = j; k < i; k++)
				{
					PrElem = PrElem + LowTr[i * MSide + k] * LInv[k * MSide + j];
				}
				LInv[i * MSide + j] = -PrElem / LowTr[i * MSide + i];
			}
		}
	}
	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			for (int k = 0; k < MSide; k++)
			{
				MatOut[i * MSide + j] = MatOut[i * MSide + j] + (Prom[i * MSide + k] * LInv[j * MSide + k]);
			}
		}
	}
	LInv.clear();
}
void CholDecomBlock(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& L, int& MTileSide, int& BlSide)
{
	std::vector<float>Prom;
	Prom.resize(BlSide * BlSide);

	int count = 0;
	for (int j = 0; j < MTileSide; j++)
	{
		for (int i = count; i < MTileSide; i++)
		{
			WithZeroVec(Prom);
			if (i == j)
			{
				if ((i != 0) && (j != 0))
				{
					CalcMat(A, L, Prom, BlSide, i, j, MTileSide);
					ChDecOneBl(Prom, L[j * MTileSide + i], BlSide);
				}
			}
			else
			{
				if (i > j)
				{
					CalcMat(A, L, Prom, BlSide, i, j, MTileSide);
					LowInvTr(L[j * MTileSide + j], Prom, L[i * MTileSide + j], BlSide);
				}
			}
		}
		count++;
	}
	Prom.clear();
}
void WithZeroMatr(std::vector<std::vector<float>>& matrix)
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[i].size(); j++)
		{
			matrix[i][j] = 0;
		}
	}
}
void CheckSol(std::vector<std::vector<float>>& L, std::vector<float>& Inp, int& MSide, int& BlSide)
{
	std::vector<float>LVec, Test;
	LVec.resize(MSide * MSide);
	Test.resize(MSide * MSide);
	CombBlock(L, LVec, MSide, BlSide);
	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			for (int k = 0; k < MSide; k++)//{}
			{
				Test[i * MSide + j] = Test[i * MSide + j] + (LVec[i * MSide + k] * LVec[j * MSide + k]);
			}
		}
	}
	if (Inp.size() == Test.size())
	{
		//for (int i = 0; i < Inp.size(); i++)
		//{
			//std::cout << Inp[i] << " " << Test[i] << std::endl;
		//}
		float TestPer1 = 0.0, SubPer = 0.0;

		for (int i = 0; i < Inp.size(); i++)
		{
			SubPer = Inp[i] - Test[i];
			TestPer1 = TestPer1 + std::pow(SubPer, 2);
		}
		TestPer1 = std::sqrt(TestPer1);
		std::cout << TestPer1 << std::endl;

		float TestPer2 = 0.0;
		for(int i = 0; i < Inp.size(); i++)
		{
			TestPer2 = TestPer2 + std::pow(Inp[i],2);
		}
		std::cout<<TestPer2<<std::endl;
		std::cout << "ForobN = "<<TestPer1 / TestPer2 << std::endl;

	}
	LVec.clear();
	Test.clear();
}
int main()
{
	std::vector<float>Inp;

	//int nTile = 16, TileSize = 4;
	//int MSide = 8, BlSide = 2;
	//int MTileSide = 4;

	int nTile = 25, TileSize = 4;
	int MSide = 10, BlSide = 2;
	int MTileSide = 5;


	//int nTile = 4, TileSize = 25;
	//int MSide = 10, BlSide = 5;
	//int MTileSide = 2;

	//int nTile = 25, TileSize = 25;
	//int MSide = 25, BlSide = 5;
	//int MTileSide = 5;


	//Заполнение исходной матрицы 
	Inp.resize(MSide * MSide);
	float PowPer = 0.0;
	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			PowPer = -((i - j) * (i - j));
			Inp[i * MSide + j] = std::exp(PowPer);
			std::cout << Inp[i * MSide + j] << " ";
		}
		std::cout << "\n";
	}

	std::vector<std::vector<float>>A(nTile, std::vector<float>(TileSize));
	std::vector<std::vector<float>>L(nTile, std::vector<float>(TileSize));

	//Разделение матрицы А на блоки
	SplitToBlock(A, Inp, MSide, BlSide);
	std::cout << "Matrix A: " << std::endl;
	PrMatrix(A);

	//Заполнение матрицы L нулями
	WithZeroMatr(L);
	std::cout << "Matrix L: " << std::endl;
	PrMatrix(L);

	//Разложение первого блока
	ChDecOneBl(A[0], L[0], BlSide);
	std::cout << "Matrix L(First block): " << std::endl;
	PrMatrix(L);

	//Разложение по блокам
	CholDecomBlock(A, L, MTileSide, BlSide);
	std::cout << "PrBlockMatrix L: " << std::endl;
	PrMatrix(L);

	//Проверка
	std::cout << "CheckSol: " << std::endl;
	CheckSol(L, Inp, MSide, BlSide);

	Inp.clear();
	//A.clear();
	//L.clear();
}
