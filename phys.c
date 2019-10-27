/* Simple 2D physics for circles using ASCII graphics
	-original NCurses code from "Game Programming in C with the Ncurses Library"
	 https://www.viget.com/articles/game-programming-in-c-with-the-ncurses-library/
	 and from "NCURSES Programming HOWTO"
	 http://tldp.org/HOWTO/NCURSES-Programming-HOWTO/
	-Physics code and algorithms from "How to Create a Custom 2D Physics
	 Engine: The Basics and Impulse Resolution"
	 https://gamedevelopment.tutsplus.com/tutorials/how-to-create-a-custom-2d-physics-engine-the-basics-and-impulse-resolution--gamedev-6331
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ncurses.h>
#include <mpi.h>     /* For MPI functions, etc */ 


// used to slow curses animation
#define DELAY 50000

// number of balls
#define POPSIZE 20

// ball radius, all circles have the same radius
#define RADIUS 1.0

// indicate if balls collide or not
#define COLLIDE 1
#define NOCOLLIDE 0

// restitution controls how bouncey the objects will be
#define RESTITUTION 0.9

// object mass
#define MASS 1.0

// maximum screen size, both height and width
#define SCREENSIZE 100

// ball location (x,y,z) and velocity (vx,vy,vz) in BallArray[][]
#define BX 0
#define BY 1
#define VX 2
#define VY 3

// ball location ERRONEOUS 1010101010101010
#define ERRONEOUS -21846

// MPI tags
# define PHYS_MPI_SEND_MASTER   20
# define PHYS_MPI_SEND_SLAVE    21
# define PHYS_MPI_LOCK_STATUS   22

// maximum screen dimensions
int MaxY = 0, MaxX = 0;

// location and velocity of ball
float BallArray[POPSIZE][4];

// change in velocity is stored for each ball (x,y,z)
float BallUpdate[POPSIZE][2];

// record the split positiosn for this job.
int* Splits;

// shared lock to break loop
int Lock;


// initialize the splitting positions
void initSplits(int commSize){

    // variables
    double splitSize;       // the ideal split size for each section
    int currentSize;        // tracks the size of the current split
    int currentSplit;       // current split position
    int processSize;        // the current size of the process


    // allocate memory for splits
    processSize = commSize - 1;
    Splits = malloc(sizeof(int) * (processSize + 1));


    // assign the first and last memebrs of the split
    Splits[0] = 0;
    Splits[processSize] = POPSIZE;
    
    // assign the rest of the split positions
    splitSize =     POPSIZE / processSize; // calculate the split size
    currentSize =   0;                  // assign the current size to 0 to start
    currentSplit =  1;                  // assign currentSplit to 1 before starting the loop
    
    
    // calculate the splits by summing approxamate splits
    for(int i = 0; i < POPSIZE; i++) {

        // increase the current size
        currentSize++;

        // currentSize 
        if(currentSize >= splitSize){
        
            // update the split index the position
            Splits[currentSplit++] = i;

            // reset the current split
            currentSize = 0;
        }
    }
    
    Splits[processSize] = POPSIZE;
}

// initialize the balls
void initBalls() {
    
    // variables
    int i;

    // calculate initial random locations for each ball, 
    // scaled based on the screen size
    for (i = 0; i < POPSIZE; i++) {

        BallArray[i][BX] = (float)(random() % SCREENSIZE);
        BallArray[i][BY] = (float)(random() % SCREENSIZE);
        BallArray[i][VX] = (float)((random() % 5) - 2);
        BallArray[i][VY] = (float)((random() % 5) - 2);
        BallUpdate[i][BX] = 0.0;
        BallUpdate[i][BY] = 0.0;
    }
}

// draw the balls
int drawBalls() {

    // variables
    int c, i;
    float multx, multy;

    // update screen maximum size
    getmaxyx(stdscr, MaxY, MaxX);


    // used to scale position of balls based on screen size
    multx = (float)MaxX / SCREENSIZE;
    multy = (float)MaxY / SCREENSIZE;


    // clear the ncurses
    clear();

    // display balls
    for (i = 0; i < POPSIZE; i++){
        
        //char* nerm;
        //nerm = malloc(sizeof(char) * 4);
        //sprintf(nerm, "%d", i);
        
        mvprintw(
            (int)(BallArray[i][BX] * multy),
            (int)(BallArray[i][BY] * multx), "o");
    }

    // refresh and delay
    refresh();
    usleep(DELAY);

    // read keyboard and exit if 'q' pressed
    c = getch();
    if (c == 'q')
        return 0;
    else
        return 1;
}

// determine if two balls in BallArray collide
// calcualte distance between the two balls and compare to the
//	sum of the radii
// use balli and ballj to identify elements in BallArray[]
int ballCollision(int balli, int ballj){

    // variables
    float distance;
    float radiiSum;


    // Pythagorean distance
    distance = sqrtf(powf((BallArray[balli][BX] - BallArray[ballj][BX]), 2) + 
        powf((BallArray[balli][BY] - BallArray[ballj][BY]), 2));
    radiiSum = RADIUS + RADIUS;
   
    // if the sum of the two radii is less than the distance
    // between the balls then they collide, otherwise they
    // do not collide
    if (distance < radiiSum)
        return (COLLIDE);
    else
        return (NOCOLLIDE);
}

// calculate the dot product between two vectors
float dotProduct(float x1, float y1, float x2, float y2) {
    
    // return the dot product
    return (x1 * x2 + y1 * y2);
}

// calculate effects of collision between BallArray[i][] and
// BallArray[j][] where i and j are the parameters to the function
void resolveCollision(int i, int j, int lowerBound, int upperBound) {

    // variables
    float rvx, rvy;
    float nx, ny;
    float distance;
    float vnormal;
    float impulse;
    float ix, iy;


    // calculate relative velocity
    rvx = BallArray[j][VX] - BallArray[i][VX];
    rvy = BallArray[j][VY] - BallArray[i][VY];

    // calculate collision normal
    nx = BallArray[j][BX] - BallArray[i][BX];
    ny = BallArray[j][BY] - BallArray[i][BY];

    // Pythagorean distance
    distance = sqrtf(powf((BallArray[j][BX] - BallArray[i][BX]), 2) + 
        powf((BallArray[j][BY] - BallArray[i][BY]), 2));

    if (distance == 0)
        return;

    nx = nx / distance;
    ny = ny / distance;

    // Calculate relative velocity in terms of the normal direction
    vnormal = dotProduct(rvx, rvy, nx, ny);

    // Do not resolve if velocities are separating
    if (vnormal > 0)
        return;

    // Calculate impulse scalar
    impulse = -(1 + RESTITUTION) * vnormal;
    impulse /= ((1 / MASS) + (1 / MASS));

    // Apply impulse
    ix = impulse * nx;
    iy = impulse * ny;

    if(lowerBound < i && i < upperBound){

        BallUpdate[i][BX] = BallArray[i][VX] - ((1 / MASS) * ix);
        BallUpdate[i][BY] = BallArray[i][VY] - ((1 / MASS) * iy);
    }

    if(lowerBound < j && j < upperBound){

        BallUpdate[j][BX] = BallArray[j][VX] + ((1 / MASS) * ix);
        BallUpdate[j][BY] = BallArray[j][VY] + ((1 / MASS) * iy);
    }
}

// calculate the movement of the balls
void moveBalls(int lowerBound, int upperBound) {
    
    // update velocity of balls based upon collisions
    // compare all balls to all other circles using two loops
    for (int i = lowerBound; i < upperBound; i++)
        for (int j = 0; j < POPSIZE; j++)
            if(i != j)
                if (ballCollision(i, j) == COLLIDE)
                    resolveCollision(i, j, lowerBound, upperBound);
    

    // move balls by calculating updating velocity and position
    for (int i = lowerBound; i < upperBound; i++) {

        // update velocity for each ball
        if (BallUpdate[i][BX] != 0.0) {
            BallArray[i][VX] = BallUpdate[i][BX];
            BallUpdate[i][BX] = 0.0;
        }
        if (BallUpdate[i][BY] != 0.0) {
            BallArray[i][VY] = BallUpdate[i][BY];
            BallUpdate[i][BY] = 0.0;
        }

        // enforce maximum velocity of 2.0 in each axis
        // done to make it easier to see collisions
        if (BallArray[i][VX] > 2.0)
            BallArray[i][VX] = 2.0;
        if (BallArray[i][VY] > 2.0)
            BallArray[i][VY] = 2.0;

        // update position for each ball
        BallArray[i][BX] += BallArray[i][VX];
        BallArray[i][BY] += BallArray[i][VY];

        // if ball moves off the screen then reverse velocity so it bounces
        // back onto the screen, and move it onto the screen
        if (BallArray[i][BX] > (SCREENSIZE - 1)) {
         
            BallArray[i][VX] *= -1.0;
            BallArray[i][BX] = SCREENSIZE - 1.5;
        }
        if (BallArray[i][BX] < 0.0) {

            BallArray[i][VX] *= -1.0;
            BallArray[i][BX] = 0.5;
        }
        if (BallArray[i][BY] > (SCREENSIZE - 1)) {

            BallArray[i][VY] *= -1.0;
            BallArray[i][BY] = SCREENSIZE - 1.5;
        }
        if (BallArray[i][BY] < 0.0) {

            BallArray[i][VY] *= -1.0;
            BallArray[i][BY] = 0.5;
        }
    }

    // mark all others before sending
    /*for(int i = 0; i < lowerBound; i++)
        for(int j = 0; j < 4; j++)
            BallArray[i][j] = ERRONEOUS;

    for(int i = upperBound; i < POPSIZE; i++)
        for(int j = 0; j < 4; j++)
            BallArray[i][j] = ERRONEOUS;*/
}

void PRIMARY_updateBallArray(float data[POPSIZE][4]){

    for(int i = 0; i < POPSIZE; i++)
        for(int j = 0; j < 4; j++)
            if(data[i][j] != ERRONEOUS)
                BallArray[i][j] = data[i][j];
}

int main(int argc, char *argv[]) {

    // MPI variables
    int commSize;       // Number of processes
    int currentRank;    // current process rank

    // variables
    int lowerBound;     // lower bound for this rank
    int upperBound;     // upper bound for this rank


    // initialize MPI
    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);       // get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);    // get current rank of this process

    // initialize upper and lower bounds
    lowerBound = -1;
    upperBound = -1;


    // initialize ncurses
    initscr();
    noecho();
    cbreak();
    timeout(0);
    curs_set(FALSE);

    // Global var `stdscr` is created by the call to `initscr()`
    getmaxyx(stdscr, MaxY, MaxX);

    // rank 0 is in charge of init
    initBalls();

    // set the lock to 1 (running)
    Lock = 1;
    
    // initialize splits each time
    initSplits(commSize);

    // initialize data for a split
    if(currentRank != 0){
        
        // assign upper and lower boud based on rank
        lowerBound = Splits[currentRank - 1];
        upperBound = Splits[currentRank];
    }


    // draw and move balls using ncurses
    // open MPI cannot sync barriers on while(true) loops for some reason, this 
    for(int i = 0; i < 200 && Lock; i++) {

        // calculate the movements on all other processes
        if(currentRank != 0){
            
            // recieve updated ball positions
            if(i != 0) // dont prefrom the first time
                MPI_Recv(
                    &BallArray, 
                    POPSIZE * 4, 
                    MPI_FLOAT, 
                    0, 
                    PHYS_MPI_SEND_SLAVE, 
                    MPI_COMM_WORLD, 
                    MPI_STATUS_IGNORE);
                    
            // calculate the movements
            moveBalls(lowerBound, upperBound);

            // send data to process 0
            MPI_Send(
                BallArray, 
                POPSIZE * 4, 
                MPI_FLOAT, 
                0, 
                PHYS_MPI_SEND_MASTER, 
                MPI_COMM_WORLD);
        }

        // wait until other processes are complete before rendering
        MPI_Barrier(MPI_COMM_WORLD);
        
        // render ncurses on process 0
        if(currentRank == 0){

            // receive data from others
            for(int j = 1; j < commSize; j++){
                
                float tempBallArray[POPSIZE][4];
                
                MPI_Recv(
                    &tempBallArray, 
                    POPSIZE * 4, 
                    MPI_FLOAT, 
                    j, 
                    PHYS_MPI_SEND_MASTER, 
                    MPI_COMM_WORLD, 
                    MPI_STATUS_IGNORE);

                PRIMARY_updateBallArray(tempBallArray);
            }
            
            // draw the balls
            Lock = drawBalls();

            // send data back
            for(int j = 1; j < commSize; j++){
            
                MPI_Send(
                    BallArray, 
                    POPSIZE * 4, 
                    MPI_FLOAT, 
                    j, 
                    PHYS_MPI_SEND_SLAVE, 
                    MPI_COMM_WORLD);

                MPI_Send(
                    &Lock, 
                    1, 
                    MPI_INT, 
                    j, 
                    PHYS_MPI_LOCK_STATUS, 
                    MPI_COMM_WORLD); 
            }
        }

        // barrier
        MPI_Barrier(MPI_COMM_WORLD);

        // receive for shutdown
        if(currentRank != 0)
            MPI_Recv(
                &Lock, 
                1,
                MPI_INT, 
                0, 
                PHYS_MPI_LOCK_STATUS, 
                MPI_COMM_WORLD, 
                MPI_STATUS_IGNORE);
    }

    endwin();

    // print out the splits at the end
    if(currentRank == 0)
        for(int i = 1; i < commSize; i++)
            printf("%d %d %d\n", i, Splits[i - 1], Splits[i]);


    // Shut down MPI
    MPI_Finalize(); 

    return 0;
}
