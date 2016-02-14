#include "CMatrixBlock.h"
#include <iostream>
#include <string.h>

using namespace std;

MatrixBlock::MatrixBlock(){

}

MatrixBlock::MatrixBlock(int nrows, int ncols)
:nrows_(nrows), ncols_(ncols)
{
	if (nrows <= 0 || ncols <= 0){
		//throw "Matrix constructor has non-posize size";
	}



	data_=new double[nrows_*ncols_];
	
	memset(data_, 0, nrows_*ncols_*sizeof(double));
}

double& MatrixBlock::operator() (int i, int j){
	if (i >= nrows_ || i < 0 || j >= ncols_ || j < 0)
		throw "Matrix subscript out of bounds";
	return data_[nrows_*i + j];
}

double MatrixBlock::operator() (int i, int j) const{
	if (i >= nrows_ || i < 0 || j >= ncols_ || j < 0)
		throw "Matrix subscript out of bounds";
	return data_[nrows_*i + j];
}


const double* MatrixBlock::get_pointer(int n) const{
	return data_+n;
}

MatrixBlock::~MatrixBlock(){
	delete[] data_;
}
