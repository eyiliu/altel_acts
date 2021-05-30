#pragma once

#include <mutex>
#include <future>

#include <cstdio>

#include "DataFrame.hh"
#include "TcpConnection.hh"

#include "myrapidjson.h"

namespace altel{
  class Layer{
  public:
    std::future<uint64_t> m_fut_async_watch;
    std::vector<DataFrameSP> m_vec_ring_ev;
    DataFrameSP m_ring_end;

    uint64_t m_size_ring{200000};
    std::atomic<uint64_t> m_count_ring_write;
    std::atomic<uint64_t> m_hot_p_read;
    uint64_t m_count_ring_read;
    bool m_is_async_watching{false};

    uint64_t m_extension{0};

    //status variable:
    std::atomic<uint64_t> m_st_n_tg_ev_now{0};
    std::atomic<uint64_t> m_st_n_ev_input_now{0};
    std::atomic<uint64_t> m_st_n_ev_output_now{0};
    std::atomic<uint64_t> m_st_n_ev_bad_now{0};
    std::atomic<uint64_t> m_st_n_ev_overflow_now{0};
    std::atomic<uint64_t> m_st_n_tg_ev_begin{0};

    uint64_t m_st_n_tg_ev_old{0};
    uint64_t m_st_n_ev_input_old{0};
    uint64_t m_st_n_ev_bad_old{0};
    uint64_t m_st_n_ev_overflow_old{0};

    std::string m_st_string;
    std::mutex m_mtx_st;

    uint32_t m_tg_expected{0};
    uint32_t m_flag_wait_first_event{true};

    bool m_isDataAccept{false};
    std::unique_ptr<TcpConnection> m_conn;
public:
    Layer(const std::string& name, const std::string& host, short int port)
      :m_name(name), m_host(host), m_port(port){
    }


    ~Layer();

    void start();
    void stop();
    void init();
    int perConnProcessRecvMesg(void* pconn, msgpack::object_handle& oh);

    DataFrameSP GetNextCachedEvent();
    DataFrameSP& Front();
    void PopFront();
    uint64_t Size();
    void ClearBuffer();

    std::string GetStatusString();
    uint64_t AsyncWatchDog();

    std::string m_name;
    std::string m_host;
    short int m_port;

  };
}
