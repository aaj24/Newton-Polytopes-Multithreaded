//
//  main.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 2/13/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "PolyDecomp.hpp"
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;

int main(int argc, const char * argv[]) {
    
    Polygon ConvexHull;
    int threads = atoi(argv[2]), level = 3;
    
    fstream file;
    file.open(argv[1]);
    if (!file.is_open())
        return 9;
        
    coord_t x, y;
    while (file >> x >> y)
            ConvexHull.push_back({x, y}); // read input
    
    file.close();

//    ConvexHull.push_back({0,0});ConvexHull.push_back({2,0});ConvexHull.push_back({2,1});ConvexHull.push_back({0,1});
    
//    ConvexHull = convex_hull(ConvexHull);
    
    double start, end;
    start = omp_get_wtime();
    
    init(ConvexHull);
    
    end = omp_get_wtime();
    
    cout << "Preprocessing time: " << end - start << endl;
    
    bool result = false;
    int u = -1;
    
    switch (threads) {
        case 0:
            info();
            break;
        case -1:
            result = PolyDecomp();
            break;
        case 1:
            result = PolyDecompDFS(level);
            break;
        default:
            result = MultiThreadedPolyDecompDFS(threads, level);
            break;
    }
    
    cout << boolalpha << result << endl;
    
//    start = omp_get_wtime();
//    threads == 1 ? r = PolyDecompNum() : r = MultiTHreadedPolyDecompNum(threads);
//    threads == 1 ? u = PolyDecompNum().u : u = MultiThreadedPolyDecompNumX(threads).u;
//    r = PolyDecompNumX();
    
//    if (threads == 1)
//        u = PolyDecompNum().u;
//    else if (threads > 1)
//        u = MultiThreadedPolyDecompNum(threads).u;
//    else if (threads < 0)
//        u = MultiThreadedPolyDecompNumX(-threads).u;
    
//    end = omp_get_wtime();
//    cout << "u = " << (u/100000) << endl;
//    cout << "Alg2 time: " << end - start << " seconds, pTp = " << (end - start)*threads << endl;
    
//    MultiThreadedSummandRecovery(2);
//    SummandRecovery();
    
//    start = omp_get_wtime();
//    if (threads == 1)
//       SummandRecovery();
//    else if (threads > 1)
//        MultiThreadedSummandRecovery(threads);
//    else if (threads < 0)
//        MultiThreadedSummandRecoveryX(-threads);
//    end = omp_get_wtime();
//    cout << "Summand Recovery time: " << end - start << " seconds" << endl;
    
    clean();
    
    return 0;
}
