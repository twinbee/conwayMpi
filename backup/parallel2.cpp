//////////////////////////////////////////////////////////////////////////////
//
// Name: Matthew Bennett
// Date: 11/6/2005
// Class: CSC730, Parallel Programming. Zhang
// Assignment: Conway's Game of Life, parallel collective version
// Description:  This is a very simple parallel implementation of Conway's
//  game of life, the 2D binary cellular automaton. 
//
// Please see the included conway.doc or conway.pdf for an explanation
//
// The workhorse is the update(petridish) function, which implments CA rules:
//  1. living cells with < 2 neighbors die at t+1
//  2. living cells with 2-3 neighbors survive at t+1 
//  3. living cells w > 3 neighbors die from overcrowding
//  4. dead cells w 3 neighbors become living
//
// This version of the code uses collective operations SCATTER and GATHER in 
//  place of the old SEND and RECEIVE.
//
// Number of processes p must satisfy p = 2^k + 1 
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// This copy of the code uses ghost_update
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//#include <cstdlib>

#include <vector>
#include <string>

using namespace std;

#include "mpi.h"

typedef vector < int >   petrirow;
typedef vector <petrirow> petridish;

typedef int adish [100][100];

petridish load(string filename);
void print(petridish lifeboard);

//work with arrays for more efficiency

void arrayify(petridish input, adish & result);

void print(adish, int sx, int sy );
void print(adish, int sx, int sy, ofstream & outfile);

void update (adish & board, int sx, int sy, int step);
void ghost_update (adish & board, int sx, int sy, int step);

void merge (adish a, adish & b, int step, int sx, int sy);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
 MPI::Init(argc, argv);

 int size = MPI::COMM_WORLD.Get_size();
 //size corresponds to the total number of processes
 
 int rank = MPI::COMM_WORLD.Get_rank();
 //rank is the process number in [0..size-1] for this process

 petridish board, localdish;
 ofstream outfile;

 int iterations = 1;
 int cur_iter = 0;

if (rank == 0)
{ 
 char choice = 'n';
 string filename;
 
 while(choice != 'y' && iterations > 0)
 { 
  cout << "The dimensions of the file determine the board's dimensions.\n";
  cout << "Please type a valid LIFE filename[simple.txt]: ";
  cin >> filename;
 
  board = load(filename);

  print(board);
  cout << "Is this the board you want? (y/n) ";
  cin.get(choice); cin.get(choice);
  cout << "How many iterations? ";  
  cin >> iterations;
 }

 outfile.open((filename + ".out").c_str());
}

 int boardh = board.size();
 MPI::COMM_WORLD.Bcast(&boardh,      1, MPI::INT, 0);
 //since P0 loaded the board, it should know the num of rows. BCAST from P0

 MPI::COMM_WORLD.Bcast(&iterations,   1, MPI::INT, 0);
 //since P0 got user input. BCAST from P0 
 
 int step = int((boardh)/(size-1));
 
 if (rank==0) while (step < 4) //assure that each process gets at least 4 rows
 //if this is not true, the algorithm loses efficiency
 // we sacrifice slave processors to achieve this result
 {
  cout << "Reducing the number of slave processes by 1\n";
  size--;
  step = int(boardh/(size-1));
  MPI::COMM_WORLD.Bcast(&size,   1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&step,   1, MPI::INT, 0); 
 }

 adish aboard;
 //convert the board into an array and also save the x size and y size
 arrayify(board, aboard);
 int aboard_sy = board.size();
 int aboard_sx = board[0].size();

//COMMUNICATION PATTERN IMPLEMENTATION FOLLOWS/////////////////////////////

 double time, start;
 
 start = MPI::Wtime(); //check and store the clock for timing purposes
 
 while ( ++cur_iter <= iterations )
 {
  if (rank==0) //master process
  {

  }
  else //slave process receive 1 or 2 ghost rows
  {   

  } 
 
 }

 time = MPI::Wtime() - start; //calculuate running time
 
// END OF IMPLEMENTATION OF COMMUNICATION PATTERN ////////////////////////////

 if (rank == 0)
 {
  cout << "Entire COLLECTIVE process took " << time << " seconds" << endl;
  outfile << "Entire COLLECTIVE process took " << time << " seconds" << endl;
  outfile.close();
 }

 MPI::Finalize();
 return 0;
}

//////////////////////////////////////////////////////////////////////////////
// loads a file into the petridish
//////////////////////////////////////////////////////////////////////////////

petridish load(string filename)
{
  petridish lifeboard;
  char buff = 'O';

  ifstream infile;

  infile.open(filename.c_str());
  if (infile.fail())
  {
   cout << "Could not open file. Using simple.txt" << endl;
   return load("simple.txt");
  }

 while (!infile.eof())
 {
  petrirow * lifeboardline = new petrirow;

  infile.get(buff);
  while (buff != '\n')
  {
   if (buff == '0')      (*lifeboardline).push_back(false);
   else if (buff == '1') (*lifeboardline).push_back(true);
   else if (buff == 'O') (*lifeboardline).push_back(false);
   else if (buff == 'X') (*lifeboardline).push_back(true);   
   infile.get(buff);
  }

  lifeboard.push_back(*lifeboardline);
 }
 
 lifeboard.pop_back();
 
 // sanity check: make sure the board contains even rows
 for (int i=1; i<lifeboard.size()-1; i++)
 {
  if (lifeboard[i].size() != lifeboard[i-1].size())
  {
   cout << "Invalid file dimensions. Board must be rectangular!"
   << " Loading simple.txt instead" << endl;
   return load("simple.txt");
  }
 }

  infile.close(); 

  return lifeboard;
}

////////////////////////////////////////////////////////////////////
// print simply prints the current state of any board passed to it
///////////////////////////////////////////////////////////////////
void print(petridish board)
{
  for (int y = 0; y < board.size(); y++)
  { for (int x = 0; x < board[y].size(); x++)
    { 
     if (board[y][x] == 0) cout << "_";
     else if (board[y][x] == 1) cout << "O";
    }
    cout <<  y << endl; 
  }
}

void print(adish board, int sx, int sy )
{
  for (int y = 0; y < sy; y++)
  { for (int x = 0; x < sx; x++)
    { 
     if (board[y][x] == 0) cout << "_";
     else if (board[y][x] == 1) cout << "O";
    }
    cout <<  y << endl; 
  }
}

////////////////////////////////////////////////////////////////////
// print simply prints the current state of any board passed to it
///////////////////////////////////////////////////////////////////
void print(adish board, int sx, int sy, ofstream & outfile)
{
  for (int y = 0; y < sy; y++)
  { for (int x = 0; x < sx; x++)
    { 
     if (board[y][x] == 0) outfile << " ";
     else if (board[y][x] == 1) outfile << "O";
    }
    outfile << endl; 
  }
}

////////////////////////////////////////////////////////////////////
// update the board according to the CA rules (simple serialized ver)
///////////////////////////////////////////////////////////////////
void update (adish & board, int sx, int sy, int step)
{
  int neighbors = 0;
  adish newboard;
  
  for (int i=0; i<sy; i++)
  for (int j=0; j<sx; j++)
   newboard[i][j] = board[i][j];
  
  for (int y = step; y < sy; y += 1)
  { for (int x = 0; x < sx; x++)
    { 
     neighbors = 0;
     if (x-1 > 0 && y-1 > 0 && x-1 < sx && y-1 < sy
          && board[y-1][x-1] == 1) neighbors++;
     if (x   > 0 && y-1 > 0 && x   < sx && y-1 < sy
          && board[y-1][x  ] == 1) neighbors++;          
     if (x+1 > 0 && y-1 > 0 && x+1 < sx && y-1 < sy
          && board[y-1][x+1] == 1) neighbors++;

     if (x-1 > 0 && y   > 0 && x-1 < sx && y   < sy
          && board[y  ][x-1] == 1) neighbors++;
     if (x+1 > 0 && y   > 0 && x+1 < sx && y   < sy
          && board[y  ][x+1] == 1) neighbors++;

     if (x-1 > 0 && y+1 > 0 && x-1 < sx && y+1 < sy
          && board[y+1][x-1] == 1) neighbors++;
     if (x   > 0 && y+1 > 0 && x   < sx && y+1 < sy
          && board[y+1][x  ] == 1) neighbors++;
     if (x+1 > 0 && y+1 > 0 && x+1 < sx && y+1 < sy
          && board[y+1][x+1] == 1) neighbors++;
     
      switch (neighbors) //implement CA rules
      {
       case 0:
       case 1:
        newboard[y][x] = 0; //die from loneliness
        break;
       case 2:
        break; //do nothing. survive or stay dead
       case 3:
        newboard[y][x] = 1; //come to life
        break;
       default: //four or more
        newboard[y][x] = 0;  //die from overcrowding
      }
    }

  }

  for (int i=0; i<sy; i++)
  for (int j=0; j<sx; j++)
   board[i][j] = newboard[i][j];
}

////////////////////////////////////////////////////////////////////
// update only ghost rows of the board
///////////////////////////////////////////////////////////////////
void ghost_update (adish & board, int sx, int sy, int step)
{
  int neighbors = 0;
  adish newboard;
  
  for (int i=0; i<sy; i++)
  for (int j=0; j<sx; j++)
   newboard[i][j] = board[i][j];
  
  for (int y = step; y < sy; y += step)
  { for (int x = 0; x < sx; x++)
    { 
     neighbors = 0;
     if (x-1 > 0 && y-1 > 0 && x-1 < sx && y-1 < sy
          && board[y-1][x-1] == 1) neighbors++;
     if (x   > 0 && y-1 > 0 && x   < sx && y-1 < sy
          && board[y-1][x  ] == 1) neighbors++;          
     if (x+1 > 0 && y-1 > 0 && x+1 < sx && y-1 < sy
          && board[y-1][x+1] == 1) neighbors++;

     if (x-1 > 0 && y   > 0 && x-1 < sx && y   < sy
          && board[y  ][x-1] == 1) neighbors++;
     if (x+1 > 0 && y   > 0 && x+1 < sx && y   < sy
          && board[y  ][x+1] == 1) neighbors++;

     if (x-1 > 0 && y+1 > 0 && x-1 < sx && y+1 < sy
          && board[y+1][x-1] == 1) neighbors++;
     if (x   > 0 && y+1 > 0 && x   < sx && y+1 < sy
          && board[y+1][x  ] == 1) neighbors++;
     if (x+1 > 0 && y+1 > 0 && x+1 < sx && y+1 < sy
          && board[y+1][x+1] == 1) neighbors++;
     
      switch (neighbors) //implement CA rules
      {
       case 0:
       case 1:
        newboard[y][x] = 0; //die from loneliness
        break;
       case 2:
        break; //do nothing. survive or stay dead
       case 3:
        newboard[y][x] = 1; //come to life
        break;
       default: //four or more
        newboard[y][x] = 0;  //die from overcrowding
      }
    }

  }

  for (int i=0; i<sy; i++)
  for (int j=0; j<sx; j++)
   board[i][j] = newboard[i][j];
}

//////////////////////////////////////////////////////////////////////////////
// merges the normal update rows with the ghost update rows
//////////////////////////////////////////////////////////////////////////////
void merge (adish a, adish & b, int sx, int sy, int step)
{
 for (int i=0; i<sy; i++)
  if (i%step)
  {
   for (int j=0; j<sx; j++) b[i][j] = a[i][j];
  }  
}

//////////////////////////////////////////////////////////////////////////////
// arrayify allows for scatter / gather to work
//////////////////////////////////////////////////////////////////////////////

void arrayify(petridish input, adish & result)
{
 for (int y = 0; y < input.size(); y++)
  for (int x = 0; x < input[y].size(); x++)  
   result[x][y] = input[x][y];
}
