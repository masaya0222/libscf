#ifndef _libiwl_iwl_hpp_
#define _libiwl_iwl_hpp_


namespace scf {

typedef short int Label;
typedef double Value;

#define IWL_INTS_PER_BUF 2980

class IWL {
    int fd_;
    int ints_per_buf_;
    int bufszc_;
    double cutoff_;
    int lastbuf_;
    int inbuf_;
    int idx_;
    Label *labels_;
    Value *values_;

public:
    IWL();
    ~IWL();

    int &fd() { return ints_per_buf_; }
    int &ints_per_buffer() { return ints_per_buf_; }
    int &buffer_size() { return bufszc_; }
    double &cutoff() { return cutoff_; }
    int &last_buffer() { return lastbuf_; }
    int &buffer_count() { return inbuf_; }
    int &index() { return idx_; }
    Label *labels() { return labels_; }
    Value *values() { return values_; }

    void init_write();
    void init_read();
    void init(const char* filename, bool for_write);

    void close_buf();
    void fetch();
    void put();

    void write(int p_first, int q_first, int r_first, int s_first, int plen, int qlen, int rlen, int slen, double* value_buf);
    void flush(int lastbuf);

};

}

#endif