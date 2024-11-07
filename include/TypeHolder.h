#ifndef TYPEHOLDER_H
#define TYPEHOLDER_H
#include <cstring>
template<typename T>
class type {
public:
    T getDefaultEps(T val){
        if(std::strcmp(typeid(T).name(),"f")==0)return 1.19E-7;
        return 2.20E-16;
    }
};
#endif