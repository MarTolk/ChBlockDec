#include<iostream>
#include<vector>
#include<cmath>
#include<lapacke.h>
#include<cblas.h>
//========================================================================================================================
void FillMat(std::vector<float>&Inp, int&MSide)
{
    Inp.resize(MSide * MSide);
    float PowPer = 0.0;
    for (int i = 0; i < MSide; i++)
    {
	for (int j = 0; j < MSide; j++)
	{
	    PowPer = -((i - j) * (i - j));
	    Inp[i * MSide + j] = std::exp(PowPer);
	}
    }
}
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
void PrVecToMatr(std::vector<float>&Vec, int&MSide)
{
    for(int i = 0; i < MSide; i++)
    {
	for(int j = 0; j < MSide; j++)
	{	
	     std::cout<<Vec[i*MSide+j]<<' ';
	}
	std::cout<<"\n";
   }
}
void WithZeroVec(std::vector<float>& Vec)
{
    for (int i = 0; i < Vec.size(); i++)
    {
	Vec[i] = 0.0;
    }
}
void WithZeroMatr(std::vector<std::vector<float>>& matrix)
{
     for (int i = 0; i < matrix.size(); i++)
     {
	for (int j = 0; j < matrix[i].size(); j++)
	{
	    matrix[i][j] = 0.0;
	}
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
void ChDecOneBl(std::vector<float>&SubA, std::vector<float>&SubL, int&BlSide)
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
					Sum = Sum + (SubL[j*BlSide+k] * SubL[j*BlSide+k]);
				}
				SubL[j*BlSide+j] = std::sqrt(SubA[j*BlSide+j] - Sum);
			}
			else 
			{
				for (int m = 0; m < j; m++)
				{
					Sum = Sum + (SubL[j*BlSide+m] * SubL[i*BlSide+m]);
				}
				SubL[i*BlSide+j] = (SubA[i*BlSide+j]-Sum);
				SubL[i*BlSide+j] = SubL[i*BlSide+j] / SubL[j*BlSide+j];
			}
		}
	}
}

void CalcMat(std::vector<std::vector<float>>&A, std::vector<std::vector<float>>&L,
	std::vector<float>& MatOut, int& BlSide, int&IndI, int&IndJ, int& MTileSide)
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
			PrElem = 0.0;
		}
	}
	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			for (int k = 0; k < MSide; k++)
			{
				MatOut[i * MSide + j] = MatOut[i * MSide + j] + (Prom[i * MSide + k] * LInv[j* MSide + k]);
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
void ChSolWithLap(std::vector<float>&Inp, std::vector<float>&LVec, std::vector<float>&Test, int&MSide)
{
    Test = Inp;
    int lda = 0, info = 0;
    lda = std::max(1,MSide);
    info = LAPACKE_spotrf(LAPACK_ROW_MAJOR, 'L', MSide, Test.data(), lda);
    if(LVec.size() == Test.size())
	{
		for(int i = 0; i < LVec.size(); i++)
		{
			if(LVec[i] == 0.0)
			{
				Test[i] = 0.0;
			}
		}
		std::cout<< "First matr "<<std::endl;
		PrVecToMatr(LVec, MSide);
		std::cout<< "LAP matr "<<std::endl;
		PrVecToMatr(Test, MSide);
		for(int i = 0; i < MSide; i++)
		{
			for(int j = 0; j < MSide; j++)
			{	
				std::cout<<LVec[i*MSide+j] <<" "<<Test[i*MSide+j]<<std::endl;
			}
		}
	}
}

void CheckSolFrobN(std::vector<float>&Inp, std::vector<float>&LVec, std::vector<float>&Test, int&MSide)
{
	WithZeroVec(Test);
	for (int i = 0; i < MSide; i++)
	{
		for (int j = 0; j < MSide; j++)
		{
			for (int k = 0; k < MSide; k++)
			{
				Test[i * MSide + j] = Test[i * MSide + j] + (LVec[i * MSide + k] * LVec[j * MSide + k]);
			}
		}
	}
	if (Inp.size() == Test.size())
	{
		float TestPer1 = 0.0, SubPer = 0.0;

		for (int i = 0; i < Inp.size(); i++)
		{
			SubPer = Inp[i] - Test[i];
			TestPer1 = TestPer1 + (SubPer*SubPer);
		}
		TestPer1 = std::sqrt(TestPer1);
		std::cout << TestPer1 << std::endl;

		float TestPer2 = 0.0;
		for(int i = 0; i < Inp.size(); i++)
		{
			TestPer2 = TestPer2 + (Inp[i]*Inp[i]);
		}
		TestPer2 = std::sqrt(TestPer2);
		std::cout<<TestPer2<<std::endl;
		std::cout << "FrobN = "<<TestPer1 / TestPer2 << std::endl;
	}
}
int main()
{
	std::vector<float>Inp;

	int nTile = 4, TileSize = 9;
	int MSide = 6, BlSide = 3;
	int MTileSide = 2;

	// int nTile = 9, TileSize = 4;
	// int MSide = 6, BlSide = 2;
	// int MTileSide = 3;

	std::vector<std::vector<float>>A(nTile, std::vector<float>(TileSize));
	std::vector<std::vector<float>>L(nTile, std::vector<float>(TileSize));

	std::vector<float>Test;
	std::vector<float>LVec;
	Test.resize(MSide * MSide);
	LVec.resize(MSide * MSide);
	
	//Заполнение исходной матрицы 
	FillMat(Inp, MSide);
	std::cout << "Fill Inp Matrix: " << std::endl;
	PrVecToMatr(Inp, MSide);

	//Разделение матрицы А на блоки
	SplitToBlock(A, Inp, MSide, BlSide);
	std::cout << "Matrix A(Block): " << std::endl;
	PrMatrix(A);
	
	//Заполнение матрицы L нулями
	WithZeroMatr(L);
	std::cout << "Matrix L(Zero): " << std::endl;
	PrMatrix(L);

	//Разложение первого блока
	ChDecOneBl(A[0], L[0], BlSide);
	std::cout << "Matrix L(First block): " << std::endl;
	PrVec(L[0]);

	//Разложение по блокам
	CholDecomBlock(A, L, MTileSide, BlSide);
	std::cout << "PrBlockMatrix L: " << std::endl;
	PrMatrix(L);

	//Проверки
	CombBlock(L, LVec, MSide, BlSide);
	
	//Проверка LAPACK для всей матрицы
	std::cout << "ChSolWithLap: " << std::endl;
	ChSolWithLap(Inp, LVec, Test, MSide);
	
	//Проверка FrobN
	std::cout << "CheckSolFrobN: " << std::endl;
	CheckSolFrobN(Inp, LVec, Test,MSide);

	Inp.clear();
	LVec.clear();
	Test.clear();
}
