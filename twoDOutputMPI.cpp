#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <mpi.h>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <iomanip>

double x(int i);                    
double y(int i);
double exactSolution(double x, double y);    
double S(double x, double y);                 
std::pair<int, int> procToLoc(int l);

double a = 0.0, b=1.0, c = 0.0, d = 1.0;
const int m=81;
const int n=81;
const int P = 9;
const int rootP = 3;
double dx=(b-a)/(m-1);
double dy=(d-c)/(n-1);
void printTotalGrid(double** grid);
double time1, time2;



double exactSolution(double x, double y) {
    return sin(x) + sin(y);
}

double x(int j, int l) {
    return (a + procToLoc(l).first * dx * m/rootP)  + j * dx;
}
double y(int k, int l) {
    return (c + procToLoc(l).second * dy * n/rootP) + k * dx;
}
double S(double x, double y) {
        return -sin(x) - sin(y);
}
enum Direction {
    UP = 0,
    DOWN,
    LEFT,
    RIGHT,
    NUMDIRS
};



double tolerance=1E-15;
int maxIterations=1000;


typedef struct _processor {
    /* We have two extra points here to cover the borders of what this processor is responsible for. */
    double** Un;
    double** Unp1;
    double** localn;
    double** localnp1;
    std::pair<int, int> loc;
    int l;
    bool commUp, commDown, commLeft, commRight;
} proc;

void printGrid2(proc me);
double** totalGrid(proc);
void procNextIter(proc* me);
std::pair<double, double> procGetErrors(proc me);
proc procInit(int l);
void procTick(proc*);
double finalGrid [m][n];
void printGrid(proc me);

int main(int argc, char* argv[])
{
    MPI_Init(NULL, NULL);

    int commSz;
    int myRank;
    int iterations = 0;
    

    MPI_Comm_size(MPI_COMM_WORLD, &commSz); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    proc me = procInit(myRank);
    
    
    // (iterationError, solutionError) 
    std::pair<double, double> localErrors;
    double globalIterError = 1., globalSolnError = 1.;
    time1 = MPI_Wtime(); 

    while (globalIterError > tolerance && iterations < maxIterations) {

        procTick(&me);

        localErrors = procGetErrors(me);
        procNextIter(&me);
        MPI_Allreduce(&localErrors.first, &globalIterError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&localErrors.second, &globalSolnError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        iterations++;
        double** grid = totalGrid(me);
     //   if (me.l == 0) printTotalGrid(grid); 

    }
    time2 = MPI_Wtime();

    if (me.l == 0) {
        // Output:
        std::cout                                                              << std::endl << std::endl;
        std::cout<< "-------------------------------------------------------"               << std::endl;
        std::cout<< "SUMMARY:"                                                 << std::endl << std::endl;
        std::cout<< "The error between two iterates is "    << globalIterError << std::endl << std::endl;
        std::cout<< "The maximum error in the solution is " << globalSolnError               << std::endl;
        std::cout<< "The number of iterations is " << iterations << std::endl;
        std::cout <<"The time taken is " << (time2 - time1) << std::endl;
        std::cout<< "-------------------------------------------------------"  << std::endl << std::endl;
    }
    /*
    std::ofstream outputFile;
    outputFile.open("output.csv");
    std::ofstream outputExact;
    outputExact.open("outputExact.csv");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            outputFile << x(i) << ',' << y(j) << ',' << Unp1[i][j] << std::endl; 
            outputExact << x(i) << ',' << y(j) << ',' << exactSolution(x(i), y(j)) << std::endl; 
        }
    }
    */
    MPI_Finalize();
    return 0;
}

// Given an 'l', return a location
std::pair<int, int> procToLoc(int l) {
    std::div_t div = std::div(l, rootP);
    return std::make_pair(div.rem, div.quot);
}
// Given a location, returns an 'l'
int locToProc(std::pair<int, int> loc) {
    return loc.first + (rootP) * loc.second;
}

proc procInit(int l) {
    proc toRet;
    toRet.Un = new double*[2 + m/rootP];
    toRet.Unp1 = new double*[2 + m/rootP];
    for (int i = 0; i < (2 + m/rootP); i++) {
        toRet.Un[i] = new double[2 + m/rootP];
        toRet.Unp1[i] = new double[2 + m/rootP];
    }
    
    
    toRet.localn = &(toRet.Un[1]);
    toRet.localnp1 = &(toRet.Unp1[1]);
   
   
    for (int i = 0; i < m/rootP + 2; i++) {
        for (int j = 0; j < n/rootP + 2; j++) {
            toRet.Un[i][j] = 0.0;
            toRet.Unp1[i][j] = 0.0;
        }
    }

     
    toRet.loc = procToLoc(l);
    toRet.commLeft = (toRet.loc.first - 1 >= 0);
    toRet.commRight = (toRet.loc.first + 1 < rootP);
    toRet.commUp = (toRet.loc.second - 1 >= 0);
    toRet.commDown = (toRet.loc.second + 1 < rootP);
    toRet.l = l;
    
    // If we are not communicating with the leftmost proc, then we set all of the values to our left 
    if (!toRet.commLeft) {
        // i actually represents i - 1 in our real grid
        for (int i = 1; i < n/rootP + 1; i++) {
            // Our grid is padded on each side by 1 row and 1 column, so we have to subtract 1 from the x value here
            toRet.Un[1][i] = exactSolution(x(0, l), y(i - 1, l));
            toRet.Unp1[1][i] = exactSolution(x(0, l), y(i - 1, l));
        }
    }
    
    if (!toRet.commRight) {
        for (int i = 1; i < n/rootP + 1; i++) {
            // Our grid is padded on each side by 1 row and 1 column, so we have to subtract 1 from the x value here
            toRet.Un[m/rootP][i] = exactSolution(x(m/rootP - 1, l), y(i - 1, l));
            toRet.Unp1[m/rootP][i] = exactSolution(x(m/rootP - 1, l), y(i - 1, l));
        }
    }
    
    if (!toRet.commUp) {
        for (int i = 1; i < m/rootP + 1; i++) {
            // Our grid is padded on each side by 1 row and 1 column, so we have to subtract 1 from the x value here
            toRet.Un[i][1] = exactSolution(x(i - 1, l), y(0, l));
            toRet.Unp1[i][1] = exactSolution(x(i - 1, l), y(0, l));
        }
    }
    
    if (!toRet.commDown) {
        for (int i = 1; i < n/rootP + 1; i++) {
            // Our grid is padded on each side by 1 row and 1 column, so we have to subtract 1 from the x value here
            toRet.Un[i][n/rootP] = exactSolution(x(i - 1, l), y(n/rootP - 1, l));
            toRet.Unp1[i][n/rootP] = exactSolution(x(i - 1, l), y(n/rootP - 1, l));
        }
    }
  
    return toRet;
}

void procTick(proc* theProc) {
    double** Unp1 = theProc->localnp1;
    double** Un = theProc->localn;
    MPI_Request outReqs[NUMDIRS];
    MPI_Request inReqs[NUMDIRS];
    double tempUp[n/rootP];
    double tempDown[n/rootP];
    double tempLeft[n/rootP];
    double tempRight[n/rootP];
    // This line only works if n = m 
    double recvBuff[NUMDIRS][n/rootP];
    
    // Non Blocking Sends and then Recieves
    if (theProc->commUp) {
        // SEND
        int dest = locToProc(std::make_pair(theProc->loc.first, theProc->loc.second - 1));
        int j = 1;
        for (int i = 0; i < m/rootP; i++) {
            tempUp[i] = Un[i][j];
        }
        MPI_Isend(tempUp, m/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &outReqs[UP]);
        // RECEIVE
        MPI_Irecv(recvBuff[UP], m/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &inReqs[UP]);
    }
    if (theProc->commDown) {
        // SEND
        int j = m/rootP;

        int dest = locToProc(std::make_pair(theProc->loc.first, theProc->loc.second + 1));

        for (int i = 0; i < m/rootP; i++) {
            tempDown[i] = Un[i][j];
        }


        MPI_Isend(tempDown, m/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &outReqs[DOWN]);
        // RECEIVE
        MPI_Irecv(recvBuff[DOWN], m/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &inReqs[DOWN]);
    }

    if (theProc->commLeft) {
        // SEND
        // Initialize our temporary left array if we need it.
        int dest = locToProc(std::make_pair(theProc->loc.first - 1, theProc->loc.second));
        int i = 0;
        for (int j = 1; j < m/rootP + 1; j++) {
            tempLeft[j - 1] = Un[i][j];
        }
        MPI_Isend(tempLeft, n/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &outReqs[LEFT]);
        // RECEIVE
        MPI_Irecv(recvBuff[LEFT], n/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &inReqs[LEFT]);
    }
    if (theProc->commRight) {
        // SEND
        int dest = locToProc(std::make_pair(theProc->loc.first + 1, theProc->loc.second));
        int i = m/rootP - 1;
        for (int j = 1; j < n/rootP + 1; j++) {
            tempRight[j - 1] = Un[i][j];
        }
        MPI_Isend(tempRight, n/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &outReqs[RIGHT]);

        // RECEIVE
        MPI_Irecv(recvBuff[RIGHT], n/rootP, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &inReqs[RIGHT]);

    }




    // Update interior points
    for (int i = 1; i < m/rootP - 1; i++) {
        for (int j = 2; j < n/rootP; j++) {
            Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
        }
    }

    
    // Update exterior points if we need to. If we don't, they're boundary conditions
    if (theProc->commUp) {
        // Block on reception 
        MPI_Wait(&inReqs[UP], MPI_STATUS_IGNORE);

        
        int j = 0;
        // Update edges
        for (int i = 0; i < m/rootP; i++) {
            Un[i][j] = recvBuff[UP][i];
        }


        // We are updating the first row MINUS THE CORNERS !!!
        j = 1;
        for (int i = 1; i < m/rootP - 1; i++) {

            // There are negative indexes here, but that's fine because our grid is padded.
            Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
        }
    }


    if (theProc->commDown) {
        // Block on reception 
        MPI_Wait(&inReqs[DOWN], MPI_STATUS_IGNORE);

        
        int j = n/rootP + 1;
        // Update edges
        for (int i = 0; i < m/rootP; i++) {
            Un[i][j] = recvBuff[DOWN][i];
        }

        // Updating the last row
        j = n/rootP;
        for (int i = 1; i < m/rootP - 1; i++) {
            // There are negative indexes here, but that's fine because our grid is padded.
            Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
        }
    }


    if (theProc->commLeft) {
        // Block on reception 
        MPI_Wait(&inReqs[LEFT], MPI_STATUS_IGNORE);

        
        int i = -1;
        // Update edges
        for (int j = 1; j < n/rootP + 1; j++) {
            Un[i][j] = recvBuff[LEFT][j - 1];
        }

        i = 0;
        for (int j = 2; j < n/rootP; j++) {
            // There are negative indexes here, but that's fine because our grid is padded.
            Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
        }
    }
    
    if (theProc->commRight) {
        // Block on reception 
        MPI_Wait(&inReqs[RIGHT], MPI_STATUS_IGNORE);

        
        // UPDATE EDGES OUTSIDE OUR BOX
        int i = m/rootP;
        for (int j = 1; j < n/rootP + 1; j++) {
            Un[i][j] = recvBuff[RIGHT][j - 1];
        }

        i = m/rootP - 1;
        for (int j = 2; j < n/rootP; j++) {
            Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
        }
    }




    // Update corners. NO need for corner updatei f they aren't communicating as they are boundary conditions.
    // Top left
    if (theProc->commUp and theProc->commLeft) {
        int i = 0, j = 1;
        Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
    }
    // Top right
    if (theProc->commUp and theProc->commRight) {
        int i = m/rootP - 1, j = 1;
        Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
    }
    // Bottom left
    if (theProc->commLeft and theProc->commDown) {
        int i = 0, j = n/rootP;
        Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
    }
    if (theProc->commRight and theProc->commDown) {
        int i = m/rootP - 1, j = n/rootP;
        Unp1[i][j] = (dy*dy*dx*dx*S(x(i, theProc->l), y(j - 1, theProc->l)) - \
                             dy*dy*(Un[i -1][j] + Un[i + 1][j]) - \
                             dx*dx*(Un[i][j -1] + Un[i][j+1])) / (-2*dy*dy - 2*dx*dx);
    }


    // Now we block for the outgoing messages so that our arrays don't go out of scope (lol)
    if (theProc->commUp)
        MPI_Wait(&outReqs[UP], MPI_STATUS_IGNORE);
    if (theProc->commDown)
        MPI_Wait(&outReqs[DOWN], MPI_STATUS_IGNORE);
    if (theProc->commLeft)
        MPI_Wait(&outReqs[LEFT], MPI_STATUS_IGNORE);
    if (theProc->commRight)
        MPI_Wait(&outReqs[RIGHT], MPI_STATUS_IGNORE);
 

}


//
// Yeah there's like a two line solution, whatever
void procNextIter(proc* me) {
    for (int i = 1; i < m/rootP + 1; i++) {
        for (int j = 1; j < m/rootP +1; j++) {
            me->Un[i][j] = me->Unp1[i][j];
        }
    }
}

std::pair<double, double> procGetErrors(proc me) {
    double iterationError = 0.0, solutionError = 0.0;
    std::pair<double, double> errorLoc;
    for(int i = 1; i< m/rootP + 1; i++){
        for (int j = 1; j < n/rootP + 1; j++) {
            double localError = fabs(me.Unp1[i][j] - me.Un[i][j]);
            double localSolutionError = fabs(me.Unp1[i][j] - exactSolution(x(i - 1, me.l), y(j - 1, me.l)) );
            if (localError > iterationError) {
                iterationError = localError;
                errorLoc = std::make_pair(x(i -1, me.l), y(j - 1, me.l));
            }

            solutionError = (localSolutionError > solutionError ? localSolutionError : solutionError);

         }
    }
    return std::make_pair(iterationError, solutionError);
}

void printGrid(proc me) {
    std::cout << "Printing grid for proc: " << me.l << std::endl;
    for (int i = 0; i < m/rootP + 2; i++) {
        for (int j = 0; j < n/rootP + 2; j++) {
            std::cout << std::setw(10) << me.Un[j][i] << " ";
        }
        std::cout << std::endl;
    }
}
void printGrid2(proc me) {
    std::cout << "Printing grid for proc: " << me.l << std::endl;
    for (int i = 0; i < m/rootP + 2; i++) {
        for (int j = 0; j < n/rootP + 2; j++) {
            std::cout << std::setw(10) << me.Unp1[j][i] << " ";
        }
        std::cout << std::endl;
    }
}

double** totalGrid(proc me) {
    double** first = new double*[m];
    for (int i = 0; i < m; i++) 
        first[i] = new double[n];

    if (me.l == 0) {
        
        for (int i = 0; i < m/rootP; i++) {
            for (int j = 0; j < m/rootP; j++) {
                first[i][j] = me.Unp1[i + 1][j + 1];
            }
        }


        for (int i = 1; i < P; i++) {
            std::pair<int, int> procLoc = procToLoc(i); 
            for (int j = 0; j < m/rootP; j++) {
                MPI_Recv(&first[m/rootP * procLoc.first + j][procLoc.second * m/rootP], m/rootP, MPI_DOUBLE, i, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            }
        }
        return first;
    }
    else {
        for (int j = 0; j < m/rootP; j++) {
            MPI_Send(&me.Unp1[j + 1][1], m/rootP, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
        }
        return NULL;
    }
}


void printTotalGrid(double** grid) {
    std::cout << "Printing grid" << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(12) << grid[j][i] << " ";
        }
        std::cout << std::endl;
    }
}

            






