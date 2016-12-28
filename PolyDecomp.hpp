//
//  PolyDecomp.hpp
//  CMPS 299
//
//  Created by Ali Jbaily on 3/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#ifndef PolyDecomp_hpp
#define PolyDecomp_hpp

#include <iostream>
#include <set>
#include <map>
#include <stack>
#include "convex_hull.hpp"
#include <math.h>
#include <algorithm>
#include <vector>
#include "omp.h"
#include "priority_update.hpp"
#include <climits>
#include "ArrayX.cpp"

typedef struct {
    int k, i, u, x, y;
} buffer_t;

typedef struct {
    unsigned int u;
    set<pair<int, int>> S;
} IntAndSet;

typedef struct {
    unsigned int u;
    ArrayX<std::pair<int, int>> a;
} IntAndX;

typedef struct {
    unsigned int u;
    map<Point, IntAndSet> *A;
} Alg2Type;

typedef struct {
    IntAndX** A;
    unsigned int u;
} Alg2TypeX;

typedef struct {
    unsigned int u;
    ArrayX<std::pair<int, int>> a;
} IntAndZ;

typedef struct {
    IntAndX** A;
    unsigned int u;
} Alg2TypeZ;

struct edge {
    int n;
    Point e;
};

struct PointAndInt{
    Point p;
    int i;
};

typedef vector<Point> Polygon;

int gcd (int, int);
void info();
bool PolyDecomp();
bool PolyDecompArray();
bool PolyDecompDFS(int level);
vector<PointAndInt> getLayer(int);
bool MultiThreadedPolyDecompDFS(int threads, int level);
Alg2TypeX PolyDecompNum();
Alg2TypeX MultiThreadedPolyDecompNum(int);
Alg2TypeX MultiThreadedPolyDecompNumX(int);
vector<stack<Point>>* SummandRecovery();
vector<stack<Point>>* SummandRecoveryX();
vector<stack<Point>>* MultiThreadedSummandRecovery(int);
vector<stack<Point>>* MultiThreadedSummandRecoveryX(int);
bool belongToIP(Point*);
void printer();
void init(Polygon);
void clean();
void initAtomic();
void deleteAtomic();
void printPoint(Point *);

#endif /* PolyDecomp_hpp */
