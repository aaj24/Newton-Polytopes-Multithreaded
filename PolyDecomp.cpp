//
//  PolyDecomp.cpp
//  CMPS 299
//
//  Created by Ali Jbaily on 3/8/16.
//  Copyright © 2016 Ali. All rights reserved.
//

#include "PolyDecomp.hpp"
#include <stdlib.h>
#include <iomanip>

Point** IP, *v0, u;
vector<Point> latticePoints;
Polygon ConvexHull;

intmax_t ymin = INT32_MAX, ymax = 0, xmin = INT32_MAX, xmax = 0;
int m, atomicLevel;
atomic<int> *I;
map<Point, atomic<int>> flags;
edge* edge_sequence;

IntAndX** newA;

atomic<bool> B;

vector<stack<Point>> *Summands;
stack<Point> *list;


void printPoint(Point *p){
    cout << p->x << ", " << p->y << endl;
}

void allocateMemory(int* size, Point* array) {
    *size <<= 1;
    Point *oldptr = array;
    array = new Point[*size];
    memcpy(array, oldptr, sizeof(Point)*(*size >> 1));
    delete [] oldptr;
}

void info(){
    int minN = INT32_MAX, maxN = 0;
    double avgN = 0, count, ratio;
    for (int i = 0; i < m; i++) {
        int n = edge_sequence[i].n;
        if (n > maxN) maxN = n;
        else if (abs(n) < minN) minN = abs(n);
        avgN += n;
    }
    
    for (coord_t y = 0; y <= ymax; y++)
    for (coord_t x = 0; x <= xmax; x++) {
        Point p = {x,y};
        if (belongToIP(&p))
        count++;
    }
    
    ratio = count/((xmax-xmin+1)*(ymax-ymin+1));
    avgN /= m;
    cout << m << "\t" << maxN << "\t" << minN << "\t" << avgN << "\t" << xmin << "\t" << xmax << "\t" << ymin << "\t" << ymax << "\t" << xmax - xmin << "\t" << ymax - ymin << "\t" << (ymax-ymin+1)*(xmax-xmin+1) << endl;
//    cout << count << "\t" << ratio << endl;
}

bool PolyDecomp() {
    
    double start = omp_get_wtime(), end;
    
    Point p;
    set<Point> oldSet, newSet;
    
    for (int i = 0; i < m - 1; i++) {
        for (int k = 1; k <= edge_sequence[i].n; k++) {
            p = *v0 + (edge_sequence[i].e * k);
            if ( belongToIP(&p))
            newSet.insert(p);
            for (auto u : oldSet) {
                p = u + edge_sequence[i].e * k;
                if (belongToIP(&p))
                newSet.insert(p);
            }
        }
        oldSet = newSet;
    }
    
    for (auto u : oldSet)
    for (int k = 0; k < edge_sequence[m-1].n; k++) {
        p = u + edge_sequence[m-1].e * k;
        if (belongToIP(&p))
        newSet.insert(p);
    }
    
    
    bool result = newSet.find(*v0) != newSet.end();
    end = omp_get_wtime();
    
    cout << "set size = " << newSet.size() << ", time = " << end - start << endl;
    return result;
}

void applyFlags(int limit) {
    if (limit > m) limit = m;
    set<Point> oldSet, newSet;
    Point p;
    oldSet.insert(*v0);
    newSet.insert(*v0);
    flags[*v0] = -1;
    
    for (int i = 0; i < limit; i++) {
        for (auto u : oldSet) {
            for (int k = 1; k <= edge_sequence[i].n; k++) {
                p = u + edge_sequence[i].e * k;
                if (belongToIP(&p)) {
                    newSet.insert(p);
                    flags[p] = -1;
                }
            }
        }
        oldSet = newSet;
    }
    atomicLevel = limit - 1;
}

void PolyDecompDFS(int i, Point *v, int stop) {
    if (B) return;
    int end;
    i < m - 1 ? end = edge_sequence[i].n : end = edge_sequence[i].n - 1;
    for (int k = 0; k <= end && !B; k++) {
        Point p = *v + edge_sequence[i].e*k;
        if (p == *v0 && k > 0) {B = true; return;}
        if (!belongToIP(&p)) return;
        else if (k > 0 && i + 1 < m)
        PolyDecompDFS(i+1, &p, m);
        else if(k == 0 && i + 1 < stop)
        PolyDecompDFS(i+1, &p, stop);
    }
}

void PolyDecompDFSAtomic(int i, Point *v, int stop) {
    if (B) return;
    int end;
    i < m - 1 ? end = edge_sequence[i].n : end = edge_sequence[i].n - 1;
    for (int k = 0; k <= end && !B; k++) {
        Point p = *v + edge_sequence[i].e*k;
        if (p == *v0 && k > 0) {B = true;return;}
        if (!belongToIP(&p)) return;
        else {
            if(flags[p] == -1) {
                priority_update(&flags[p], i);
                if (i + 1 <= atomicLevel)
                    PolyDecompDFSAtomic(i+1, &p, m);
                else if (i + 1 < m)
                    PolyDecompDFS(i+1, &p, m);
            }
            else if (flags[p] > i){
                
                if(k > 0) {
                    if (flags[p] < stop) {
                        stop = flags[p];
                        atomic_compare_exchange_strong(&flags[p], &stop, i);
                    }
                    if (i + 1 <= atomicLevel)
                        PolyDecompDFSAtomic(i+1, &p, stop);
                    else if (i + 1 < m)
                        PolyDecompDFSAtomic(i+1, &p, stop);
                }
                
                else if(k == 0 && i + 1 < stop) {
                    if (i + 1 <= atomicLevel)
                        PolyDecompDFSAtomic(i+1, &p, stop);
                    else if (i + 1 < m)
                        PolyDecompDFS(i+1, &p, stop);
                }
            }
            else if (k == 0 && i + 1 < stop) { // k == 0 so it can continue progressing from that point and not fail when i > flag[p]
                if (i + 1 <= atomicLevel)
                    PolyDecompDFSAtomic(i+1, &p, stop);
                else if (i + 1 < m)
                    PolyDecompDFS(i+1, &p, stop);
            }
        }
    }
}

bool PolyDecompDFS(int level) {
    double start, end;
    start = omp_get_wtime();
    applyFlags(level);
    end = omp_get_wtime();
    cout << "Applying flags: " << end - start << endl;
    start = omp_get_wtime();
    PolyDecompDFSAtomic(0, v0, m);
    end = omp_get_wtime();
    cout << "Serial DFS time: " << end - start << endl;
    return B;
}

vector<PointAndInt> getLayer(int threads) {
    
    vector<PointAndInt> oldSet, newSet;
    Point p;
    map<Point, int> mp;
    oldSet.push_back({*v0, 0});
    newSet.push_back({*v0, 0});
    threads--;
    
    
    for (int i = 0; i < m - 1; i++) {
        for (auto u : oldSet) {
            int step;
            if (edge_sequence[i].n == 1 || threads >= (edge_sequence[i].n - 1)) step = 1;
            else step = floor((edge_sequence[i].n - 1)/threads);
            for (int k = step; k <= edge_sequence[i].n; k+=step) {
                p = u.p + edge_sequence[i].e * k;
                if (belongToIP(&p) && mp[{p.x,p.y}] == 0) {
                    flags[p] = i;
                    newSet.push_back({p, i + 1});
                    threads--;
                    if (threads == 0) return newSet;
                }
            }
        }
        
        oldSet = newSet;
    }
    
    return newSet;
}

bool MultiThreadedPolyDecompDFS(int threads, int level) {
    double start, end, totTime;
    int tid;
    start = omp_get_wtime();
    applyFlags(level);
    vector<PointAndInt> points = getLayer(threads);
    end = omp_get_wtime();
    cout << "Applying flags: " << end - start << endl;
    omp_set_num_threads(threads);
    start = omp_get_wtime();
    Point p;
#pragma omp parallel shared(totTime) private(tid, p)
    {
        tid = omp_get_thread_num();
        double tstart, tend;
        tstart = omp_get_wtime();
        
//        int start = tid*(edge_sequence[0].n/threads), end;
        
//        p = *v0 + edge_sequence[0].e * start;
        
//        priority_update(&flags[p], 0);
//        PolyDecompDFSAtomic(0, &p, m);
        
        PolyDecompDFSAtomic(points[tid].i, &points[tid].p, m);
        PolyDecompDFSAtomic(0, v0, m);
        
        tend = omp_get_wtime();
        totTime += (tend - tstart);
#pragma omp critical
        {
            cout << "thread " << tid << ", at point (" << points[tid].p.x << ", " << points[tid].p.y << ") finished in " << tend - tstart << " seconds" << endl;
        }
    }
    end = omp_get_wtime();
    cout << "Multithreaded DFS time: " << end - start << ", ∑time = " << totTime << ", p×Tp = " << threads*(end - start) << totTime << endl;
    return B;
}

void getLatticePoints() {
    for (coord_t y = ymin; y <= ymax; y++)
        for (coord_t x = xmin; x <= xmax; x++) {
            Point v = {x, y};
            if (belongToIP(&v)){
                latticePoints.push_back(v);
            }
        }
}

Alg2TypeX PolyDecompNum() {
    
//    unsigned int consumed = 0;
    
    newA = new IntAndX*[ymax+1];
    IntAndX** oldA = new IntAndX*[ymax+1];
    for (int y = 0; y <= ymax; y++) {
        newA[y] = new IntAndX[xmax+1];
        oldA[y] = new IntAndX[xmax+1];
    }
    
    newA[(int)v0->y][(int)v0->x].u = 1;
    oldA[(int)v0->y][(int)v0->x].u = 1;
    
//    consumed = 2 * (ymax + 1) * (xmax + 2) * sizeof(IntAndX);
    
    for (int i = 0; i < m; i++) {
        for (int y = ymin; y <= ymax; y++)
            for (int x = xmin; x <= xmax; x++) {
                Point v = {x,y};
                if (oldA[y][x].u > 0) {
                    for (int k = 1; k <= edge_sequence[i].n; k++) {
                        Point vp = v + edge_sequence[i].e*k;
                        if (belongToIP(&vp)) {
                            int xp = (int)vp.x, yp = (int)vp.y;
                            newA[yp][xp].u += oldA[y][x].u;
                            newA[yp][xp].a.add({k,i});
//                            consumed += sizeof(int);
                        }
                    }
                }
            }
        
        for (int y = (int)ymin; y <= ymax; y++)
            for (int x = (int)xmin; x <= xmax; x++)
                oldA[y][x].u = newA[y][x].u;

//        consumed += (ymax - ymin + 1) * (xmax - xmin + 1) * sizeof(int);
        
//        cout << i << '\n';
    }
    
//    for (int y = (int)ymin; y <= ymax; y++)
//        for (int x = (int)xmin; x <= xmax; x++)
//            consumed += newA[y][x].a.consumed;
//    cout << consumed << " bytes consumed" << endl;
    
    for (int i = 0; i <= ymax; i++)
        delete [] oldA[i];
    delete[] oldA;
    
    
    cout << newA[(int)v0->y][(int)v0->x].a.size << endl;
    
    return {newA, newA[(int)v0->y][(int)v0->x].u};
}

bool between(Point *p, Point *p1, Point *p2) {
    if(p1->y < p->y && p->y < p2->y) // between
        return true;
    else if((p1->x <= p->x && p1->y == p->y) || (p->x <= p2->x && p->y == p2->y)) // extremities
        return true;
    else if(p1->x <= p->x && p->x <= p2->x && p1->y == p2->y) // on the same line
        return true;
    else
        return false;
}

int getSector(Point *p, int threads){
    for (int t = 0; t < threads; t++) {
        int start = t*((xmax-xmin)/threads) + xmin, end;
        t == threads - 1 ? end = xmax : end = ((t+1)*((xmax-xmin)/threads)) - 1 + xmin;
        if (p->x >= start && p->x <= end)
            return t;
    }
    cout << "-1\n";
    return -1;
}

Alg2TypeX MultiThreadedPolyDecompNum(int threads) {
    
    omp_set_num_threads(threads);
    
//    Matrix<pair<int, int>>::defaultSize = threads;
    
    newA = new IntAndX*[ymax+1];
    IntAndX** oldA = new IntAndX*[ymax+1];
//    mutex** mx = new mutex*[ymax+1];
    
#pragma omp parallel
    {
        int tid = omp_get_thread_num(), counter = 0;
        
#pragma omp for
        for (int y = 0; y <= ymax; y++) {
            newA[y] = new IntAndX[xmax+1];
            oldA[y] = new IntAndX[xmax+1];
        }
        
#pragma omp master
        {
            newA[(int)v0->y][(int)v0->x].u = 1;
            oldA[(int)v0->y][(int)v0->x].u = 1;
        }
        
        int start = tid*((xmax-xmin)/threads) + xmin, end;
        tid == threads - 1 ? end = xmax : end = ((tid+1)*((xmax-xmin)/threads)) - 1 + xmin;
        if (end < start) end = start;

        
        for (int i = 0; i < m; i++) {
            for (int y = ymin; y <= ymax; y++)
                for (int x = start; x <= end; x++) {
                    Point v = {x,y};
                    if (oldA[y][x].u > 0) {
                        for (int k = 1; k <= edge_sequence[i].n; k++) {
                            Point vp = v + edge_sequence[i].e*k;
                            if (belongToIP(&vp)) {
                                int xp = (int)vp.x, yp = (int)vp.y;
#pragma omp atomic
                                newA[yp][xp].u += oldA[y][x].u;
                                
//                                mx[yp][xp].lock();
                                newA[yp][xp].a.atomicAdd({k,i});
//                                newA[yp][xp].a.push_back({k, i});
//                                mx[yp][xp].unlock();
                            }
                        }
                    }
            }
            
#pragma omp barrier

            for (int y = ymin; y <= ymax; y++)
                for (int x = start; x <= end; x++)
                    oldA[y][x].u = newA[y][x].u;
#pragma omp barrier

        }
        
#pragma omp for
        for (int i = 0; i <= ymax; i++)
            delete [] oldA[i];
    }
    
    delete[] oldA;
    
    
    cout << newA[(int)v0->y][(int)v0->x].a.size << endl;
    
    return {newA, newA[(int)v0->y][(int)v0->x].u};
}

Alg2TypeX MultiThreadedPolyDecompNumX(int threads) {
    
    if (threads >= 23)
        threads = 1;
    
    omp_set_num_threads(threads);
    
    newA = new IntAndX*[ymax+1];
    IntAndX** oldA = new IntAndX*[ymax+1];
    
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        
#pragma omp for
        for (int y = 0; y <= ymax; y++) {
            newA[y] = new IntAndX[xmax+1];
            oldA[y] = new IntAndX[xmax+1];
        }
        
#pragma omp master
        {
            newA[(int)v0->y][(int)v0->x].u = 1;
            oldA[(int)v0->y][(int)v0->x].u = 1;
        }
        
        int start = tid*((xmax-xmin)/threads) + xmin, end;
        tid == threads - 1 ? end = xmax : end = ((tid+1)*((xmax-xmin)/threads)) - 1 + xmin;
        if (end < start) end = start;
        
        
        for (int i = 0; i < m; i++) {
            for (int y = ymin; y <= ymax; y++)
                for (int x = start; x <= end; x++) {
                    Point v = {x,y};
                    if (belongToIP(&v)) {
                        for (int k = 1; k <= edge_sequence[i].n; k++) {
                            Point vp = v - edge_sequence[i].e*k;
                            int xp = (int)vp.x, yp = (int)vp.y;
                            if (belongToIP(&vp) && oldA[yp][xp].u > 0) {
                                newA[y][x].u += oldA[yp][xp].u;
                                newA[y][x].a.add({k,i});
                            }
                        }
                    }
                }
            
#pragma omp barrier
            
            for (int y = ymin; y <= ymax; y++)
                for (int x = start; x <= end; x++)
                    oldA[y][x].u = newA[y][x].u;
#pragma omp barrier
            
        }
        
#pragma omp for
        for (int i = 0; i <= ymax; i++)
            delete [] oldA[i];
    }
    
    delete[] oldA;
    
    
    cout << newA[(int)v0->y][(int)v0->x].a.size << endl;
    
    return {newA, newA[(int)v0->y][(int)v0->x].u};
}

void SummandRecovery(Point p, int prv){
    int x = (int)p.x, y = (int)p.y;
    for (int j = 0; j < newA[y][x].a.size; j++) {
        int k = newA[y][x].a[j].first, i = newA[y][x].a[j].second;
        if (i < prv) {
            Point v = p - edge_sequence[i].e * k;
            list->push(v);
            if (v == *v0 && i < m - 1) Summands->push_back(*list);
            else SummandRecovery(v, i);
            list->pop();
        }
    }
}

vector<stack<Point>>* SummandRecovery(){
    list = new stack<Point>();
    Summands = new vector<stack<Point>>();
    
    double start = omp_get_wtime(), end;
    
    SummandRecovery(*v0, m+1);
    
    end = omp_get_wtime();
    
    cout << "time: " << end - start << ", " << Summands->size() << " collected" << endl;
    delete list;
    return Summands;
}

void MultiThreadedSummandRecovery(Point p, int prv, int tid){
    int x = (int)p.x, y = (int)p.y;
    for (int j = 0; j < newA[y][x].a.size; j++) {
        int k = newA[y][x].a[j].first, i = newA[y][x].a[j].second;
        if (i < prv) {
            Point v = p - edge_sequence[i].e * k;
            list[tid].push(v);
            if (v == *v0 && i < m - 1)
                Summands[tid].push_back(list[tid]);
            else MultiThreadedSummandRecovery(v, i, tid);
            list[tid].pop();
        }
    }
}

vector<stack<Point>>* MultiThreadedSummandRecovery(int threads){

    cout << "threads: " << threads << endl;
    
    double start = omp_get_wtime(), end;
    
    unsigned int size = newA[(int)v0->y][(int)v0->x].a.size;
    
    pair<int, int> *pairs = new pair<int, int>[size];
    atomic<int> *locks = new atomic<int>[size];
    list = new stack<Point>[threads];
    Summands = new vector<stack<Point>>[threads];
    
    int pos = 0;
    for(int i = 0; i < newA[(int)v0->y][(int)v0->x].a.size; i++)
        pairs[pos++] = newA[(int)v0->y][(int)v0->x].a[i];
    
    for (int i = 0; i < size; i++)
        locks[i] = -1;
    
    omp_set_num_threads(threads);
    
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for (int p = 0; p < size; p++) {
            if (locks[p] == -1)
                priority_update(&locks[p], tid);
            if (locks[p] == tid) {
                pair<int, int > pair = pairs[p];
                int k = pair.first, i = pair.second;
                Point v = *v0 - edge_sequence[i].e * k;
                list[tid].push(v);
                if (v == *v0 && i < m - 1) Summands[tid].push_back(list[tid]);
                else MultiThreadedSummandRecovery(v, i, tid);
                list[tid].pop();
            }
        }
    }
    
    end = omp_get_wtime();
    
    int totSize = 0;
    for (int t = 0; t < threads; t++)
        totSize += Summands[t].size();
    
    cout << "time: " << end - start << ", " << totSize << " collected" << endl;
    return Summands;
}

vector<stack<Point>>* MultiThreadedSummandRecoveryX(int threads){
    
    cout << "threads = " << threads << endl;
    
    double start = omp_get_wtime(), end;
    
    unsigned int size = newA[(int)v0->y][(int)v0->x].a.size;
    
    pair<int, int> *pairs = new pair<int, int>[size];
    atomic<int> *locks = new atomic<int>[size];
    list = new stack<Point>[threads];
    Summands = new vector<stack<Point>>[threads];
    
    int pos = 0;
    for(int i = 0; i < newA[(int)v0->y][(int)v0->x].a.size; i++)
        pairs[pos++] = newA[(int)v0->y][(int)v0->x].a[i];
    
    for (int i = 0; i < size; i++)
        locks[i] = -1;
    
    omp_set_num_threads(threads);
    
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp for
        for (int p = 0; p < size; p++) {
            pair<int, int > pair = pairs[p];
            int k = pair.first, i = pair.second;
            Point v = *v0 - edge_sequence[i].e * k;
            list[tid].push(v);
            if (v == *v0 && i < m - 1) Summands[tid].push_back(list[tid]);
            else MultiThreadedSummandRecovery(v, i, tid);
            list[tid].pop();
        }
    }
    
    
    end = omp_get_wtime();
    
    int totSize = 0;
    for (int t = 0; t < threads; t++)
        totSize += Summands[t].size();
    
    cout << "time: " << end - start << ", " << totSize << " collected" << endl;
    return Summands;
}

Point** CalcIP() {
    int convexSize = (int)ConvexHull.size();
    ymax = ymin = ConvexHull[0].y;
    for (int i = 1; i < convexSize; i++) {
        if (ConvexHull[i].y > ymax)
        ymax = ConvexHull[i].y;
        if (ConvexHull[i].y < ymin)
        ymin = ConvexHull[i].y;
        if (ConvexHull[i].x > xmax)
        xmax = ConvexHull[i].x;
        if (ConvexHull[i].x < xmin)
        xmin = ConvexHull[i].x;
    }
    
    IP = new Point*[ymax + 1];
    for (int i = 0; i <= ymax; i++)
    IP[i] = new Point[2];
    
    double slope[convexSize];
    for (int i = 0; i < convexSize - 1; i++)
    slope[i]=(ConvexHull[i].y-ConvexHull[i+1].y)/(ConvexHull[i].x - ConvexHull[i+1].x);
    
    for (intmax_t y = ymin; y <= ymax; y++) {
        int count = 0;
        double x[2] = {-1, -1};
        for (int p = 0; p < convexSize && count < 2; p++) {
            const int y1 = ConvexHull[p].y;
            const int y2 = ConvexHull[p+1].y;
            if (y >= min(y1, y2) && y <= max(y1, y2)) {
                if (abs(slope[p]) < 1e-2) {
                    x[0] = ConvexHull[p].x;
                    x[1] = ConvexHull[(p+1)%convexSize].x;
                    count = 2;
                }
                else {
                    double value = ((y-ConvexHull[p].y)/slope[p])+ConvexHull[p].x;
                    if (x[0] != value)
                    x[count++] = value;
                }
            }
        }
        const int x1 = (const int)x[0], x2 = (const int)x[1];
        IP[y][0].x = ceil(min(x1, x2));
        IP[y][1].x = floor(max(x1, x2));
        IP[y][0].y = IP[y][1].y = y;
    }
    
    if (IP[ymin][0].x == -1)
    IP[ymin][0].x = IP[ymin][1].x;
    if (IP[ymax][0].x == -1)
    IP[ymax][0].x = IP[ymax][1].x;
    
    return IP;
}

bool belongToIP(Point* p) {
    return p->y <= ymax && p->y >= ymin && p->x >= IP[(int)p->y][0].x && p->x <= IP[(int)p->y][1].x;
}

int gcd (int a, int b) {
    int c;
    while (a != 0) {
        c = a; a = b%a;  b = c;
    }
    return b;
}

void printer() {
    for (int i = 0; i < ymax; i++) {
        cout << "(" << IP[i][0].x << ", " << IP[i][0].y << ")" << "\t";
        cout << "(" << IP[i][1].x << ", " << IP[i][1].y << ")" << endl;
    }
}

void init(Polygon p) {
    ConvexHull = p;
    m = (int)ConvexHull.size() - 1;
    edge_sequence = new edge[m];
    
    for (int i = 0; i < m; i++) {
        coord_t dx = ConvexHull[i+1].x - ConvexHull[i].x, dy = ConvexHull[i+1].y - ConvexHull[i].y;
        int n = gcd(dx, dy);
        n > 0 ? edge_sequence[i] = {n, dx/n, dy/n} : edge_sequence[i] = {-n, dx/-n, dy/-n};
    }
    
    v0 = &ConvexHull[0];
    CalcIP();
    cout << "m = " << m << ", n0 = " << edge_sequence[0].n << endl;
}

void clean() {
    delete [] edge_sequence;
    for (int i = 0; i <= ymax; i++)
    delete [] IP[i];
    delete [] IP;
    delete newA;
}
