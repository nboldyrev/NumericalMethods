#include "linearAlgebra.h"
#include <iostream>


namespace Lab1 {//* NOTE::  файлы 
    const std::string inDir = "../tests/L1/";
    const std::string outDir = "../results/L1/";

    const std::string fileIn1 = inDir+"test1.txt";
    const std::string fileIn2 = inDir+"test2.txt";
    const std::string fileIn3 = inDir+"test3.txt";
    const std::string fileIn4 = inDir+"test4.txt";
    const std::string fileIn5 = inDir+"test5.txt";
    const std::string dopFileIn1 = inDir+"dop_test1.txt";
    const std::string dopFileIn2 = inDir+"dop_test2.txt";

    const std::string fileOut1 = outDir+"result1.txt";
    const std::string fileOut2 = outDir+"result2.txt";
    const std::string fileOut3 = outDir+"result3.txt";
    const std::string fileOut4 = outDir+"result4.txt";
    const std::string fileOut5 = outDir+"result5.txt";
    const std::string dopFileOut1 = outDir+"result6.txt";
    const std::string dopFileOut2=outDir+"result7.txt";
}

namespace Lab2 {//* NOTE:: файлы
    const std::string inDir = "../tests/L2/";
    const std::string outDir = "../results/L2/";

    const std::string fileIn1 = inDir+"test1.txt";
    const std::string fileIn2 = inDir+"test2.txt";
    const std::string fileIn3 = inDir+"test3.txt";

    const std::string fileOut1 = outDir+"result1.txt";
    const std::string fileOut2 = outDir+"result2.txt";
    const std::string fileOut3 = outDir+"result3.txt";
   
}




int main() {
    GaussSolver solG;
    SimpleIterSolver solS(0,0.0001);
    auto problem0 = solG.readLSE(Lab1::fileIn2);
    auto problem1 = solS.readLSE(Lab2::fileIn1);
    std::cout<<solG.solve(problem0)<<"\n\n";
    std::cout<<solS.solve(problem1)<<"\n\n";
    return 0;
}
