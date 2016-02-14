#ifndef __CCONN_H__
#define __CCONN_H__

#include <string.h>

class ConnList{
public:
	ConnList(int nconns) :
        nconns_(nconns)
    {
        l_ = new int[nconns_];
        r_ = new int[nconns_];
        trans_ = new double[nconns_];
        memset(l_,0,nconns_*sizeof(int));
        memset(r_,0,nconns_*sizeof(int));
        memset(trans_,0,nconns_*sizeof(double));
    }
    const int* get_l() const{ return l_; }
    const int* get_r() const { return r_; }
    const double* get_trans() const { return trans_; }
    bool InsertConnection(int n,int l,int r,double trans){
        l_[n] = l;
        r_[n] = r;
        trans_[n] = trans;
        return true;
    }
    bool InsertConnection(int* l,int* r,double* trans){
        for(int i = 0;i<nconns_;i++){
            l_[i] = l[i];
            r_[i] = r[i];
            trans_[i] = trans[i];
        }
        return true;
    }
	~ConnList(){
        delete []l_;
        delete []r_;
        delete []trans_;
    }
private:
    int nconns_;
    int* l_;
    int* r_;
    double* trans_;
};

#endif
