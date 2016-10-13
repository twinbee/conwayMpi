//////////////////////////////////////////////////////////////////////////////
//
// Name: Matthew Bennett
// Date: 11/6/2005
// Class: CSC730, Parallel Programming. Zhang
// Assignment: Conway's Game of Life, parallel BLOCKING version
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
// This version of the code uses BLOCKING operations SEND and RECV instead of
//  collective operations SCATTER and GATHER
//
// Number of processes p must satisfy p = 2^k + 1 
//
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// This copy of the code uses ghost_update
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#include <cstdlib>

#include <vector>
#include <string>

using namespace std;

#include "mpi.h"

typedef vector < int >   petrirow;
typedef vector <petrirow> petridish;

petridish load(string filename);
void print(petridish lifeboard);
void print(petridish lifeboard, ofstream & outfile);

void sendrow (const petridish & board,
                const int &k, const int &dest, const int &tag);
petrirow recvrow (const int &src, const int &tag);


petridish update (petridish board);
petridish ghost_update (petridish board, int step);

void merge (petridish a, petridish & b, int step);

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

 int iterations = 50;
 int cur_iter = 0;

if (rank == 0)
{ 
 char choice = 'n';
 string filename = "simple.txt";
 
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

//COMMUNICATION PATTERN IMPLEMENTATION FOLLOWS/////////////////////////////

 double time, start;
 
 start = MPI::Wtime(); //check and store the clock for timing purposes
 
 while ( ++cur_iter <= iterations )
 {
  // first send and receive rows (delegate from P0)

  if (rank==0) //master process
  {
   for (int i=0*step; i< 1*step+1; i++) //send p1 one ghost row at end
    sendrow(board, i, 1, cur_iter);

   for (int dest=2; dest < size-1; dest++)
   for (int i=(dest-1)*step-1; i< dest*step+1; i++)
    sendrow(board, i, dest, cur_iter); //send intermediate procs 2 ghost rows

   for (int i=(size-2)*step-1; i< (size-1)*step; i++)
    sendrow(board, i, size-1, cur_iter);  //send p-last one ghost row beginning
  }
  else //slave process receive 1 or 2 ghost rows
  {   
   localdish.clear();

   if (rank==1 || rank==size-1)
   for (int i=0; i<step+1; i++)
    localdish.push_back(recvrow(0, cur_iter));
   
   else for (int i=0; i<step+2; i++)
    localdish.push_back(recvrow(0, cur_iter));
  
  } 

  // --- MASTER AND SLAVES ARE IN SYNC DUE TO BLOCKING SEND --- 
  //accumulate the data at P0
  if (rank == 0) // (master process)
  {
   localdish.clear();
   
   for (int i=0; i <step; i++) //one ghost row
    localdish.push_back(recvrow(1, 1000+cur_iter));   
    recvrow(1, 1000+cur_iter); //throw away row     

   for (int src=2; src < size-1; src++) //two ghost rows
   {
    recvrow(src, 1000+cur_iter); //throw away row     
    for (int i=0; i<step; i++)
     localdish.push_back(recvrow(src, 1000+cur_iter));
//    recvrow(src, 1000+cur_iter); //throw away row          
   }

   recvrow(size-1, 1000+cur_iter); //throw away row
   for (int i=0; i<step; i++) //one ghost row
    localdish.push_back(recvrow(size-1, 1000+cur_iter));

   board = ghost_update(board, step);
   merge(localdish, board, step);

   cout << "--BOARD-----------Iteration: " << cur_iter << endl;
   print(board);
   print(board, outfile);
  }
  //run update at each process except for the master and return the data
  else // (slave process)
  {
   localdish = update(localdish);

   for (int i=0; i<localdish.size(); i++)
    sendrow(localdish, i, 0, 1000+cur_iter);
  }

  // --- MASTER AND SLAVES ARE IN SYNC DUE TO BLOCKING SEND --- 
  
 }

 time = MPI::Wtime() - start; //calculuate running time
 
// END OF IMPLEMENTATION OF COMMUNICATION PATTERN ////////////////////////////

 if (rank == 0)
 {
  cout << "Entire BLOCKING process took " << time << " seconds" << endl;
  outfile << "Entire BLOCKING process took " << time << " seconds" << endl;
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

////////////////////////////////////////////////////////////////////
// update only ghost rows of the board
///////////////////////////////////////////////////////////////////
petridish ghost_update (petridish board, int step)
{
  int neighbors = 0;
  petridish newboard = board; //make a copy only for sizing
  
  for (int y = step; y < board.size(); y += step)
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
// performs a blocking send on row k of board
///////////////////////////////////////////////////////////////////////////////
void sendrow (const petridish & board,
                const int &k, const int &dest, const int &tag)
{
 int array[100];
 int size = board[k].size();

 MPI::COMM_WORLD.Send(&size, 1, MPI::INT, dest, tag);
 
 for (int i=0; i<board[k].size(); i++)
  array[i] = board[k][i];

 MPI::COMM_WORLD.Send(array, size, MPI::INT, dest, tag);
}

///////////////////////////////////////////////////////////////////////////////
// performs a blocking recv on row k of board
///////////////////////////////////////////////////////////////////////////////
petrirow recvrow (const int &src, const int &tag)
{
 int array[100];
 petrirow board;
 
 int size;
 
 MPI::COMM_WORLD.Recv(&size,    1, MPI::INT, src, tag);
 
 MPI::COMM_WORLD.Recv(array, size, MPI::INT, src, tag);

 for (int i=0; i<size; i++)
  board.push_back(array[i]);
 
 return board;
}

//////////////////////////////////////////////////////////////////////////////
// merges the normal update rows with the ghost update rows
//////////////////////////////////////////////////////////////////////////////
void merge (petridish a, petridish & b, int step)
{
 if (a.size() != b.size())
 {
  cout << "ERROR, CANNOT MERGE UNLIKE VECTORS OF SIZE " 
    << a.size() << " & " << b.size() << endl;
 }
 
 for (int i=0; i<a.size(); i+=1)
 if (i%step) b[i] = a[i];
 //else b[i] = a[i];
}

