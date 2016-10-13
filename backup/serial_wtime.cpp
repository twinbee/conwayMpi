//////////////////////////////////////////////////////////////////////////////
//
// Name: Matthew Bennett
// Date: 11/6/2005
// Class: CSC730, Parallel Programming. Zhang
// Assignment: Conway's Game of Life, serialized version
// Description:  This is a very simple serialized implementation of Conway's
//  game of life, the 2D binary cellular automaton. Nothing tricky in the code.
//  The serial version is for comparison purposes with the parallelized version
//  of the code, which is written in MPI CH. 
//
// The workhorse is the update(petridish) function, which implments CA rules:
//  1. living cells with < 2 neighbors die at t+1
//  2. living cells with 2-3 neighbors survive at t+1 
//  3. living cells w > 3 neighbors die from overcrowding
//  4. dead cells w 3 neighbors become living
//
// THIS VERSION INCLUDES MPI JUST FOR BENCHMARKING. It is still serial
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#include <cstdlib>

#include <vector>
#include <string>

#include "mpi.h"

using namespace std;

typedef vector <vector < bool > > petridish;

void delay( float n );

petridish load(string filename);
void print(petridish lifeboard);
void print(petridish lifeboard, ofstream & outfile);

petridish update (petridish board);

int main(int argc, char *argv[])
{
 MPI::Init(argc, argv);

 double start = MPI::Wtime();

 string filename("simple.txt");
 petridish board;
 ofstream outfile;

 int iterations = 50;
 int cur_iter = 0;
  
 char choice = 'n';

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

 while (++cur_iter <= iterations)
 {
  board = update(board);
  for (int i=1; i<board.size(); i++) cout << endl;
  print(board);
  print(board, outfile);
  outfile << "-- Iteration number " << cur_iter <<  "  ----------------------\n";
  delay(0.2); //delay for 1 second. Remove in parallel version.
 }
 
 double total_time = MPI::Wtime() - start;
 

cout << "Entire SERIAL process took " << total_time << " seconds" << endl;
outfile << "Entire SERIAL process took " << total_time << " seconds" << endl;

 outfile.close();

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
  vector <bool> * lifeboardline = new vector <bool>;

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
     if (board[y][x] == 0) cout << " ";
     else if (board[y][x] == 1) cout << "O";
    }
    cout << endl; 
  }
}

////////////////////////////////////////////////////////////////////
// print simply prints the current state of any board passed to it
///////////////////////////////////////////////////////////////////
void print(petridish board, ofstream & outfile)
{
  for (int y = 0; y < board.size(); y++)
  { for (int x = 0; x < board[y].size(); x++)
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
petridish update (petridish board)
{
  int neighbors = 0;
  petridish newboard = board; //make a copy only for sizing
  
  for (int y = 0; y < board.size(); y++)
  { for (int x = 0; x < board[y].size(); x++)
    { 
     neighbors = 0;
     if (x-1 > 0 && y-1 > 0 && x-1 < board[y-1].size() && y-1 < board.size()
          && board[y-1][x-1] == 1) neighbors++;
     if (x   > 0 && y-1 > 0 && x   < board[y-1].size() && y-1 < board.size()
          && board[y-1][x  ] == 1) neighbors++;          
     if (x+1 > 0 && y-1 > 0 && x+1 < board[y-1].size() && y-1 < board.size()
          && board[y-1][x+1] == 1) neighbors++;

     if (x-1 > 0 && y   > 0 && x-1 < board[y].size() && y   < board.size()
          && board[y  ][x-1] == 1) neighbors++;
     if (x+1 > 0 && y   > 0 && x+1 < board[y].size() && y   < board.size()
          && board[y  ][x+1] == 1) neighbors++;

     if (x-1 > 0 && y+1 > 0 && x-1 < board[y+1].size() && y+1 < board.size()
          && board[y+1][x-1] == 1) neighbors++;
     if (x   > 0 && y+1 > 0 && x   < board[y+1].size() && y+1 < board.size()
          && board[y+1][x  ] == 1) neighbors++;
     if (x+1 > 0 && y+1 > 0 && x+1 < board[y+1].size() && y+1 < board.size()
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

 return newboard;
}


///////////////////////////////////////////////////////////////////////////////
// crude and dirty time delay function
//////////////////////////////////////////////////////////////////////////////

void delay( float n )
{
 //Sleep(1000 * n); //only works in windows / DOS
 //sleep(1000 * n); //only works in *nix 
}
