#ifndef __CMATRIXBLOCK_H__
#define __CMATRIXBLOCK_H__



using namespace std;

class MatrixBlock{
public:
	MatrixBlock();
	MatrixBlock(int nrows, int ncols);
	double& operator() (int i, int j);
	double operator() (int i, int j) const;
	bool set_x(int x);
	bool set_y(int y);
	const double* get_pointer(int n) const;
	~MatrixBlock();

private:
	int nrows_, ncols_;
	double* data_;
};


#endif