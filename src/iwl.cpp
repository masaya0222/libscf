#include "src/iwl.hpp"

#include <iostream>
#include <unistd.h>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define SYSTEM_READ ::read
#define SYSTEM_WRITE ::write

namespace scf {

IWL::IWL() {
    ints_per_buf_ = IWL_INTS_PER_BUF;
    cutoff_ = 1.e-14;
    bufszc_ = 2 * sizeof(int) + ints_per_buf_ * 4 * sizeof(Label) + ints_per_buf_ * sizeof(Value);
    lastbuf_ = 0;
    inbuf_ = 0;
    idx_ = 0;
}

void IWL::init_write() { init("test.txt", true); }
void IWL::init_read() { init("test.txt", false); }

void IWL::init(const char* filename, bool for_write) {
    ints_per_buf_ = IWL_INTS_PER_BUF;
    cutoff_ = 1.e-14;
    bufszc_ = 2 * sizeof(int) + ints_per_buf_ * 4 * sizeof(Label) + ints_per_buf_ * sizeof(Value);
    lastbuf_ = 0;
    inbuf_ = 0;
    idx_ = 0;

    labels_ = new Label[4 * ints_per_buf_];
    values_ = new Value[ints_per_buf_];

    if (for_write)
        fd_ = open(filename, O_CREAT|O_WRONLY|O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH );
    else
        fd_ = open(filename, O_RDONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH );

}

IWL::~IWL() { close_buf(); std::cout << "iwl destructor was called" << std::endl; }

void IWL::close_buf() {
    close(fd_);
    if (labels_) delete[](labels_);
    if (values_) delete[](values_);
    labels_ = nullptr;
    values_ = nullptr;
}

void IWL::fetch() {
    SYSTEM_READ( fd_, (char *)&(lastbuf_), sizeof(int));
    SYSTEM_READ( fd_, (char*)&(inbuf_), sizeof(int));
    SYSTEM_READ( fd_, (char*)labels_, inbuf_ * 4 * sizeof(Label));
    SYSTEM_READ( fd_, (char*)values_, inbuf_ * sizeof(Value));
    idx_ = 0;
}

void IWL::put() {
    SYSTEM_WRITE( fd_, (char*)&(lastbuf_), sizeof(int));
    SYSTEM_WRITE( fd_, (char*)&(inbuf_), sizeof(int));
    SYSTEM_WRITE( fd_, (char*)labels_, inbuf_ * 4 * sizeof(Label));
    SYSTEM_WRITE( fd_, (char*)values_, inbuf_ * sizeof(Value));
}

void IWL::write(int p_first, int q_first, int r_first, int s_first, int plen, int qlen, int rlen, int slen, double* value_buf) {
    int p, q, r, s, idx;
    double value;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    for(int f1 = 0, f1234=0; f1!=plen; ++f1) {
        p = f1 + p_first;
        for(int f2 = 0; f2!=qlen; ++f2) {
            q = f2 + q_first;
            for(int f3 = 0; f3 != rlen; ++f3) {
                r = f3 + r_first;
                for(int f4 = 0; f4 != slen; ++f4, ++f1234) {
                    s = f4 + s_first;
                    
                    if (std::fabs(value_buf[f1234]) < cutoff_)
                        continue;

                    idx = 4 * idx_;
                    lblptr[idx++] = (Label)p;
                    lblptr[idx++] = (Label)q;
                    lblptr[idx++] = (Label)r;
                    lblptr[idx++] = (Label)s;
                    valptr[idx_] = (Value)value_buf[f1234];
                    idx_++;
    
                    if (idx_ == ints_per_buf_) {
                        lastbuf_ = 0;
                        inbuf_ = idx_;
                        put();
                        idx_ = 0;
                    }
                }
            }
        }
    }
}

void IWL::flush(int lastbuf) {
    int idx;
    Label *lblptr;
    Value *valptr;

    inbuf_ = idx_;
    lblptr = labels_;
    valptr = values_;

    idx = 4 * idx_;

    while (idx < ints_per_buf_)
    {
        lblptr[idx++] = 0;
        lblptr[idx++] = 0;
        lblptr[idx++] = 0;
        lblptr[idx++] = 0;
        valptr[idx_] = 0.0;
        idx_++;
    }

    if (lastbuf)
        lastbuf_ = 1;
    else
        lastbuf_ = 0;

    put();
    idx_ = 0;
}

}