#pragma once
#include <iostream>
#include <memory>
#include <map>
#include <thread>
#include <mutex>

template <class c1, class c2>
struct mutex_map
{
    void insert(std::pair<c1, c2> pair);
    void insert(c1 p1, c2 p2);
    
    std::mutex mtx;
    std::map<c1, c2> map;
};

typedef mutex_map<uint16_t, uint32_t> mutex_map_16_32;
typedef mutex_map<uint16_t, uint64_t> mutex_map_16_64;

struct mutex_map_16_64_32
{
    mutex_map_16_64 imsi_map;
    mutex_map_16_32 m_tmsi_map;
}; 
typedef mutex_map_16_64_32 rnti_map;

template <class c1, class c2>
void mutex_map<c1, c2>::insert(std::pair<c1, c2> pair)
{
    std::lock_guard<std::mutex> lg(mtx);
    map.insert(pair);
}

template <class c1, class c2>
void mutex_map<c1, c2>::insert(c1 p1, c2 p2)
{
    insert(std::pair<c1, c2>(p1, p2));
}
