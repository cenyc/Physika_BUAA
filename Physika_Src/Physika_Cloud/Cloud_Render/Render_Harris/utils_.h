//
// Created by czn on 1/20/19.
//
//#pragma once

#include <string>
#include <unistd.h>
#ifndef PHYSIKA_UTILS_H
#define PHYSIKA_UTILS_H


std::string get_cwd(){
    char cwnd[100];
    return std::string(getcwd(cwnd, sizeof(cwnd)));
}


#endif PHYSIKA_UTILS_H
