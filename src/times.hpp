//
//  times.hpp
//  cDenemeler
//
//  Created by Mehmet  Nuri on 5.05.2021.
//  Copyright © 2021 Mehmet  Nuri. All rights reserved.
//

#ifndef times_hpp
#define times_hpp

#include <stdio.h>
#include <chrono>
#include <ctime>
#include <iostream>

void mainTimes();


struct Timer{
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<float> duration;
    const char* functionName;

    Timer(const char* name){
        start = std::chrono::high_resolution_clock::now();
        functionName = name;
    }

    ~Timer(){
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;

        double time_s = duration.count();
        printf("%s: Timer result: %f s \n", functionName,time_s);
    }
};


#endif /* times_hpp */
